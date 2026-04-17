"""Report engine – SolverResult를 다양한 형식으로 렌더링.

Provides three rendering modes:
  1. ``render_markdown``   – GitHub-flavoured Markdown with LaTeX math blocks
  2. ``render_exam_style`` – Exam-answer-sheet style (Unicode box-drawing characters)
  3. ``save_markdown``     – Convenience: write Markdown to a file

Unicode math helpers are included for inline equation display in plain-text
environments where a full LaTeX renderer is not available.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .base import SolverResult

# ---------------------------------------------------------------------------
# Unicode math symbol table (subset used in ME551 vibration problems)
# ---------------------------------------------------------------------------

MATH_SYMBOLS: dict[str, str] = {
    # Greek letters – lower case
    "alpha":   "α",  "beta":  "β",  "gamma":  "γ",  "delta": "δ",
    "epsilon": "ε",  "zeta":  "ζ",  "eta":    "η",  "theta": "θ",
    "lambda":  "λ",  "mu":    "μ",  "nu":     "ν",  "xi":    "ξ",
    "pi":      "π",  "rho":   "ρ",  "sigma":  "σ",  "tau":   "τ",
    "phi":     "φ",  "chi":   "χ",  "psi":    "ψ",  "omega": "ω",
    # Greek letters – upper case
    "Gamma":   "Γ",  "Delta":  "Δ",  "Theta":  "Θ",  "Lambda": "Λ",
    "Xi":      "Ξ",  "Pi":     "Π",  "Sigma":  "Σ",  "Phi":    "Φ",
    "Psi":     "Ψ",  "Omega":  "Ω",
    # Operators & special
    "sqrt":    "√",  "inf":   "∞",  "approx": "≈",  "neq":    "≠",
    "leq":     "≤",  "geq":   "≥",  "pm":     "±",  "times":  "×",
    "dot":     "·",  "deg":   "°",  "sum":    "Σ",  "int":    "∫",
    "partial": "∂",  "nabla": "∇",
    # Subscript/superscript digits (convenience)
    "_0": "₀", "_1": "₁", "_2": "₂", "_3": "₃", "_4": "₄",
    "_5": "₅", "_6": "₆", "_7": "₇", "_8": "₈", "_9": "₉",
    "^2": "²", "^3": "³",
    # Arrows
    "to":      "→",  "implies": "⟹",  "iff": "⟺",
    # Box-drawing (used internally by render_exam_style)
    "_box_heavy": "═",  "_box_light": "─",
    "_box_tri":   "▶",  "_box_check": "✓",
    "_box_warn":  "⚠",
}


# Keys whose values are large numerical arrays — omit from printed output
# (still available in the SolverResult object for plotting / export)
_SKIP_IN_DISPLAY: set[str] = {
    "response_t", "response_x", "response_q", "response_qdot",
    "transition_matrix_values",
}


def unicode_math(text: str) -> str:
    """Replace LaTeX-style escape sequences with Unicode math characters.

    Supported patterns:
        \\omega  → ω
        \\Omega  → Ω
        \\theta  → θ
        \\zeta   → ζ
        \\sqrt   → √
        (etc. – see MATH_SYMBOLS table)

    Also converts common shorthand:
        omega_n → ωₙ,  zeta  → ζ,  theta_eq → θ_eq
    """
    # Replace \\name escapes
    for name, sym in MATH_SYMBOLS.items():
        text = text.replace(f"\\{name}", sym)
    # Common shorthand replacements
    text = text.replace("omega_n", "ωₙ")
    text = text.replace("omega_d", "ωd")
    text = text.replace("omega",   "ω")
    text = text.replace("Omega",   "Ω")
    text = text.replace("theta",   "θ")
    text = text.replace("sigma",   "σ")
    text = text.replace("lambda",  "λ")
    text = text.replace("zeta",    "ζ")
    text = text.replace("alpha",   "α")
    text = text.replace("beta",    "β")
    text = text.replace("sqrt",    "√")
    text = text.replace("''",      "″")
    text = text.replace("'",       "′")
    return text


def _fmt_value(val: Any, unicode_mode: bool = True) -> str:
    """Format a dict value for display (round floats, apply unicode math)."""
    if isinstance(val, float):
        formatted = f"{val:.6g}"
    elif isinstance(val, list):
        formatted = str(val)
    else:
        formatted = str(val)
    return unicode_math(formatted) if unicode_mode else formatted


# ---------------------------------------------------------------------------
# Exam-style box constants
# ---------------------------------------------------------------------------

_BOX_WIDTH    = 55
_HEAVY_RULE   = "═" * _BOX_WIDTH
_LIGHT_RULE   = "─" * _BOX_WIDTH
_SECTION_MARK = "▶"


def _box_header(title: str) -> str:
    """Return a centred title inside heavy-rule box lines."""
    padded = f" {title} "
    side   = max(0, (_BOX_WIDTH - len(padded)) // 2)
    line   = "═" * side + padded + "═" * (_BOX_WIDTH - side - len(padded))
    return f"{_HEAVY_RULE}\n{line}\n{_HEAVY_RULE}"


def _indent(text: str, n: int = 2) -> str:
    """Indent every line of *text* by *n* spaces."""
    pad = " " * n
    return "\n".join(pad + line for line in text.splitlines())


# ---------------------------------------------------------------------------
# ReportEngine
# ---------------------------------------------------------------------------


class ReportEngine:
    """SolverResult를 읽기 좋은 여러 형식으로 변환한다.

    All methods are static so the class can be used without instantiation.

    Usage::

        result = solver.solve(params)
        print(ReportEngine.render_exam_style(result))
        ReportEngine.save_markdown(result, "report.md")
    """

    # ------------------------------------------------------------------
    # 1. Markdown renderer
    # ------------------------------------------------------------------

    @staticmethod
    def render_markdown(result: SolverResult) -> str:
        """SolverResult → GitHub-flavoured Markdown string.

        Equations wrapped in ``$$...$$`` blocks for LaTeX rendering.
        """
        lines: list[str] = []

        # Title
        lines.append(f"# {result.problem_type}")
        lines.append("")

        # Given
        lines.append("## Given")
        for key, value in result.given.items():
            lines.append(f"- **{key}**: {value}")
        lines.append("")

        # Derivation Steps
        lines.append("## Derivation Steps")
        for i, (desc, expr) in enumerate(result.steps, 1):
            lines.append(f"### Step {i}: {desc}")
            lines.append("")
            expr_str = str(expr)
            # Decide whether to wrap in math block or code block
            if any(c in expr_str for c in ("=", "→", "∑", "∫", "∂")):
                lines.append(f"$$\n{expr_str}\n$$")
            else:
                lines.append(f"```\n{expr_str}\n```")
            lines.append("")

        # Final Answer
        lines.append("## Final Answer")
        for key, value in result.final_answer.items():
            if key in _SKIP_IN_DISPLAY:
                continue
            lines.append(f"- **{key}**: `{value}`")
        lines.append("")

        # Sanity Check
        if result.sanity_check:
            lines.append("## Sanity Check")
            lines.append(result.sanity_check)
            lines.append("")

        return "\n".join(lines)

    # ------------------------------------------------------------------
    # 2. Exam-style renderer
    # ------------------------------------------------------------------

    @staticmethod
    def render_exam_style(result: SolverResult) -> str:
        """Format SolverResult as an exam answer sheet.

        Output format::

            ═══════════════════════════════════════════════════════
            ════════════  Problem: FRF Decomposition  ════════════
            ═══════════════════════════════════════════════════════

            ▶ Given:
              ...

            ▶ Solution:
              Step 1: ...
              Step 2: ...

            ▶ Final Answer:
              ...

            ▶ Verification:
              ...
            ═══════════════════════════════════════════════════════
        """
        blocks: list[str] = []

        # Header box
        blocks.append(_box_header(f"Problem: {result.problem_type}"))
        blocks.append("")

        # ── Given ──────────────────────────────────────────────
        blocks.append(f"{_SECTION_MARK} Given:")
        if result.given:
            for key, value in result.given.items():
                blocks.append(f"  {key} = {_fmt_value(value)}")
        else:
            blocks.append("  (no given parameters)")
        blocks.append("")

        # ── Solution ────────────────────────────────────────────
        blocks.append(f"{_SECTION_MARK} Solution:")
        if result.steps:
            for i, (desc, expr) in enumerate(result.steps, 1):
                header = f"  Step {i}: {unicode_math(desc)}"
                blocks.append(header)
                # Indent the expression content
                expr_lines = str(expr).splitlines()
                for line in expr_lines:
                    blocks.append(f"    {unicode_math(line)}")
                blocks.append("")
        else:
            blocks.append("  (no steps recorded)")
            blocks.append("")

        # ── Final Answer ────────────────────────────────────────
        blocks.append(f"{_SECTION_MARK} Final Answer:")
        if result.final_answer:
            for key, value in result.final_answer.items():
                if key in _SKIP_IN_DISPLAY:
                    continue
                blocks.append(
                    f"  {unicode_math(key)} = {_fmt_value(value)}"
                )
        else:
            blocks.append("  (no final answer recorded)")
        blocks.append("")

        # ── Verification ────────────────────────────────────────
        if result.sanity_check:
            blocks.append(f"{_SECTION_MARK} Verification:")
            for line in result.sanity_check.splitlines():
                blocks.append(f"  {unicode_math(line)}")
            blocks.append("")

        # Footer
        blocks.append(_HEAVY_RULE)

        return "\n".join(blocks)

    # ------------------------------------------------------------------
    # 3. Compact summary renderer
    # ------------------------------------------------------------------

    @staticmethod
    def render_summary(result: SolverResult) -> str:
        """One-page compact summary — problem type, key answers, stability verdict.

        Useful for quick console inspection without the full derivation.
        """
        lines: list[str] = []
        lines.append(_LIGHT_RULE)
        lines.append(f"  SUMMARY: {result.problem_type}")
        lines.append(_LIGHT_RULE)

        # Key final answers only
        fa = result.final_answer
        if fa:
            lines.append("  Key results:")
            for key, val in fa.items():
                lines.append(f"    {unicode_math(key):<30s} {_fmt_value(val)}")
        else:
            lines.append("  (no results)")

        if result.sanity_check:
            lines.append(_LIGHT_RULE)
            lines.append("  Checks:")
            for line in result.sanity_check.splitlines():
                lines.append(f"    {unicode_math(line)}")

        lines.append(_LIGHT_RULE)
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # 4. Exam-answer mode
    # ------------------------------------------------------------------

    @staticmethod
    def render_exam_answer(result: SolverResult) -> str:
        """SolverResult -> 서술형 답안 스타일 출력.

        기존 render_exam_style이 '계산 과정 나열'이라면,
        이 모드는 '시험 답안처럼 결론 + 근거'를 간결하게 출력.
        """
        lines: list[str] = []
        fa = result.final_answer
        pt = result.problem_type

        lines.append(f"[{pt}]")
        lines.append("")

        # --- T/F 판별 ---
        if "ConceptDB" in pt:
            verdict = fa.get("adjusted_verdict", fa.get("verdict", ""))
            rule_id = fa.get("rule_id", "")
            reason = fa.get("reason", "")
            note = fa.get("verdict_note", "")
            confidence = fa.get("confidence", "")

            lines.append(f"판정: {verdict}")
            if note:
                lines.append(f"근거: {note}")
            if reason:
                # Truncate long reasons
                reason_short = reason[:200] + "..." if len(reason) > 200 else reason
                lines.append(f"설명: {reason_short}")
            if fa.get("common_trap"):
                lines.append(f"주의: {fa['common_trap']}")
            if fa.get("conflicts"):
                lines.append("감지된 함정:")
                for c in fa["conflicts"]:
                    lines.append(f"  - {c['detail']}")

        # --- FRF ---
        elif "FRF" in pt:
            # Key subsystem info
            for key in ("subsystems", "poles", "dc_gain"):
                if key in fa:
                    lines.append(f"{key}: {_fmt_value(fa[key])}")
            if "characteristic_values" in fa:
                lines.append("특성값:")
                for cv_key, cv_val in fa["characteristic_values"].items():
                    lines.append(f"  {cv_key}: {_fmt_value(cv_val)}")
            if "forced_response" in fa:
                lines.append(f"강제응답: {_fmt_value(fa['forced_response'])}")

        # --- Lagrange ---
        elif "Lagrange" in pt:
            if "eom" in fa:
                eom = fa["eom"]
                if isinstance(eom, list):
                    for i, e in enumerate(eom):
                        lines.append(f"EOM_{i+1}: {e} = 0")
                else:
                    lines.append(f"EOM: {eom} = 0")
            if "equilibrium_points" in fa:
                lines.append(f"평형점: {fa['equilibrium_points']}")
            if "stability" in fa:
                for sr in fa["stability"]:
                    eq = sr.get("equilibrium", "")
                    stab = sr.get("stability", "")
                    lines.append(f"안정성 ({eq}): {stab}")
            if "linearized" in fa:
                for lin in fa["linearized"]:
                    lines.append(f"선형화 M={lin.get('M','')}, K={lin.get('K','')}")

        # --- Modal ---
        elif "Modal" in pt:
            if "natural_frequencies" in fa:
                lines.append(f"고유진동수: {_fmt_value(fa['natural_frequencies'])}")
            if "mode_shapes" in fa:
                lines.append(f"모드형상: {_fmt_value(fa['mode_shapes'])}")
            if "modal_response" in fa:
                lines.append(f"모달 응답: {_fmt_value(fa['modal_response'])}")

        # --- Damping ---
        elif "Damping" in pt or "damping" in pt.lower():
            if "proportional" in fa:
                lines.append(f"비례감쇠: {fa['proportional']}")
            if "damping_ratios" in fa:
                lines.append(f"감쇠비: {_fmt_value(fa['damping_ratios'])}")
            if "alpha" in fa and "beta" in fa:
                lines.append(f"alpha={fa['alpha']}, beta={fa['beta']}")

        # --- Extended ---
        elif "Extended" in pt:
            if "system_type" in fa:
                lines.append(f"시스템 유형: {fa['system_type']}")
            if "stability" in fa:
                lines.append(f"안정성: {fa['stability']}")
            if "stability_details" in fa:
                lines.append(f"상세: {fa['stability_details']}")
            if "flutter_mode_count" in fa:
                lines.append(f"Flutter 모드: {fa['flutter_mode_count']}")
            if "divergence_mode_count" in fa:
                lines.append(f"Divergence 모드: {fa['divergence_mode_count']}")

        # --- Generic fallback ---
        else:
            for key, val in fa.items():
                lines.append(f"{key}: {_fmt_value(val)}")

        # Sanity check
        if result.sanity_check:
            lines.append("")
            lines.append(f"검증: {result.sanity_check}")

        return "\n".join(lines)

    # ------------------------------------------------------------------
    # 5. File I/O helpers
    # ------------------------------------------------------------------

    @staticmethod
    def save_markdown(result: SolverResult, filepath: str | Path) -> None:
        """Save SolverResult as a Markdown file.

        Creates parent directories as needed.
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        content = ReportEngine.render_markdown(result)
        filepath.write_text(content, encoding="utf-8")

    @staticmethod
    def save_exam_style(result: SolverResult, filepath: str | Path) -> None:
        """Save SolverResult as an exam-style plain-text file (UTF-8)."""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        content = ReportEngine.render_exam_style(result)
        filepath.write_text(content, encoding="utf-8")
