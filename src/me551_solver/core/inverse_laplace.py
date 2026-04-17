"""Inverse Laplace Transform 솔버.

F(s) = N(s)/D(s) 형태의 입력을 받아:
  1. 극점(poles) 분석
  2. 부분분수 분해 (partial fraction decomposition)
  3. 각 항의 역 라플라스 변환
  4. 최종 f(t) 도출

입력 방식:
  a) 분자/분모 계수 리스트  (예: num=[1], den=[1,4,9])
  b) sympy 표현식 문자열    (예: "1/(s**2 + 4*s + 9)")
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from .base import BaseSolver, SolverResult

# Shared symbols
s = sp.Symbol("s")
t = sp.Symbol("t", positive=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _coeffs_to_poly(coeffs: list, var=s) -> sp.Poly:
    """Convert coefficient list [high ... low] to sympy Poly."""
    return sp.Poly(coeffs, var)


def _poly_to_expr(coeffs: list, var=s) -> sp.Expr:
    """Convert coefficient list to sympy expression."""
    n = len(coeffs) - 1
    expr = sp.Integer(0)
    for i, c in enumerate(coeffs):
        c_rational = sp.Rational(c).limit_denominator(10000) if isinstance(c, float) else sp.Integer(c)
        expr += c_rational * var ** (n - i)
    return expr


def _format_expr(expr: sp.Expr) -> str:
    """Pretty-format a sympy expression."""
    return str(expr).replace("**", "^").replace("*", "·")


def _format_expr_clean(expr: sp.Expr) -> str:
    """Format expression with sqrt display."""
    text = sp.pretty(expr, use_unicode=True)
    return text


def _classify_pole(pole: sp.Expr) -> str:
    """Classify a pole for display."""
    if pole.is_real:
        if pole == 0:
            return "원점 (s=0)"
        elif pole < 0:
            return "실수 음극 (stable)"
        else:
            return "실수 양극 (unstable)"
    else:
        re_part = sp.re(pole)
        if re_part < 0:
            return "복소 음극 (stable oscillatory)"
        elif re_part == 0:
            return "순허수극 (marginally stable)"
        else:
            return "복소 양극 (unstable oscillatory)"


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

class InverseLaplaceSolver(BaseSolver):
    """F(s) → f(t) 역 라플라스 변환 솔버.

    Supports:
      - Coefficient input: num=[...], den=[...]
      - Expression string: "1/(s**2 + 4*s + 9)"
      - Direct sympy expression
    """

    def solve(self, params: dict) -> SolverResult:
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # --- Parse input ---
        F_s = self._parse_input(params)

        # Display F(s)
        steps.append((
            "Given F(s)",
            f"F(s) = {_format_expr(F_s)}"
        ))
        final_answer["F_s"] = str(F_s)

        # --- Separate into numerator / denominator ---
        numer, denom = sp.fraction(sp.together(F_s))
        numer = sp.expand(numer)
        denom = sp.expand(denom)

        steps.append((
            "Numerator / Denominator",
            f"N(s) = {_format_expr(numer)}\n"
            f"D(s) = {_format_expr(denom)}"
        ))

        # --- Factor denominator & find poles ---
        denom_factored = sp.factor(denom)
        poles_raw = sp.roots(denom, s)  # {pole: multiplicity}

        pole_strs = []
        for pole, mult in sorted(poles_raw.items(), key=lambda x: (sp.re(x[0]), sp.im(x[0]))):
            classification = _classify_pole(pole)
            mult_str = f" (중복도 {mult})" if mult > 1 else ""
            pole_strs.append(f"  s = {pole}{mult_str}  — {classification}")

        steps.append((
            "Denominator factorization & poles",
            f"D(s) = {_format_expr(denom_factored)}\n\n"
            f"Poles:\n" + "\n".join(pole_strs)
        ))
        final_answer["poles"] = {str(p): int(m) for p, m in poles_raw.items()}

        # --- Partial fraction decomposition ---
        pfd = sp.apart(F_s, s)
        pfd_expanded = sp.Add.make_args(pfd)

        pfd_strs = []
        for i, term in enumerate(pfd_expanded, 1):
            pfd_strs.append(f"  Term {i}: {_format_expr(term)}")

        steps.append((
            "Partial fraction decomposition",
            f"F(s) = {_format_expr(pfd)}\n\n" +
            "\n".join(pfd_strs)
        ))
        final_answer["partial_fractions"] = str(pfd)

        # --- Inverse Laplace of each term ---
        ilt_strs = []
        f_t_total = sp.Integer(0)

        for i, term in enumerate(pfd_expanded, 1):
            f_t_term = sp.inverse_laplace_transform(term, s, t)
            # Remove Heaviside(t) wrapper if present (we assume t > 0)
            f_t_term = f_t_term.rewrite(sp.exp)
            f_t_term = f_t_term.replace(sp.Heaviside(t), sp.Integer(1))
            f_t_term = f_t_term.replace(sp.Heaviside, lambda *args: sp.Integer(1))
            f_t_term = sp.simplify(f_t_term)

            ilt_strs.append(
                f"  L⁻¹{{{_format_expr(term)}}} = {_format_expr(f_t_term)}"
            )
            f_t_total += f_t_term

        f_t_total = sp.simplify(f_t_total)
        # Try to express in terms of sin/cos instead of exp
        f_t_trig = sp.simplify(sp.trigsimp(f_t_total.rewrite(sp.cos)))

        steps.append((
            "Inverse Laplace transform (term by term)",
            "\n".join(ilt_strs)
        ))

        # --- Final result ---
        # Choose the more compact form
        f_t_display = f_t_trig if len(str(f_t_trig)) <= len(str(f_t_total)) else f_t_total

        steps.append((
            "Final result: f(t)",
            f"f(t) = {_format_expr(f_t_display)}    (t ≥ 0)"
        ))
        final_answer["f_t"] = str(f_t_display)

        # --- Also try expanded form for readability ---
        f_t_expanded = sp.expand_trig(f_t_display)
        f_t_expanded = sp.expand(f_t_expanded)
        if str(f_t_expanded) != str(f_t_display):
            steps.append((
                "Expanded form",
                f"f(t) = {_format_expr(f_t_expanded)}"
            ))
            final_answer["f_t_expanded"] = str(f_t_expanded)

        # --- Verification: initial & final value theorems ---
        sanity_parts = []
        try:
            # Initial value: f(0+) = lim_{s→∞} s·F(s)
            iv_laplace = sp.limit(s * F_s, s, sp.oo)
            iv_direct = f_t_display.subs(t, 0)
            iv_match = sp.simplify(iv_laplace - iv_direct) == 0
            sanity_parts.append(
                f"Initial value theorem: lim(s→∞) s·F(s) = {iv_laplace}, "
                f"f(0) = {iv_direct}: {'PASS' if iv_match else 'CHECK'}"
            )
        except Exception:
            pass

        try:
            # Final value: f(∞) = lim_{s→0} s·F(s)  (if all poles in LHP or at origin)
            all_stable = all(sp.re(p) <= 0 for p in poles_raw)
            if all_stable:
                fv_laplace = sp.limit(s * F_s, s, 0)
                fv_direct = sp.limit(f_t_display, t, sp.oo)
                fv_match = sp.simplify(fv_laplace - fv_direct) == 0
                sanity_parts.append(
                    f"Final value theorem: lim(s→0) s·F(s) = {fv_laplace}, "
                    f"f(∞) = {fv_direct}: {'PASS' if fv_match else 'CHECK'}"
                )
        except Exception:
            pass

        return SolverResult(
            problem_type="Inverse Laplace Transform",
            given=self._build_given(params),
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts) if sanity_parts else "",
        )

    def _parse_input(self, params: dict) -> sp.Expr:
        """Parse input into a sympy expression F(s)."""
        # Mode 1: expression string
        expr_str = params.get("expression")
        if expr_str:
            local_dict = {"s": s, "t": t}
            return sp.sympify(expr_str, locals=local_dict)

        # Mode 2: coefficient lists
        num_coeffs = params.get("numerator", params.get("num"))
        den_coeffs = params.get("denominator", params.get("den"))
        if num_coeffs is not None and den_coeffs is not None:
            numer = _poly_to_expr(num_coeffs)
            denom = _poly_to_expr(den_coeffs)
            return numer / denom

        raise ValueError(
            "Input required: either 'expression' string or "
            "'numerator'/'denominator' coefficient lists."
        )

    def _build_given(self, params: dict) -> dict:
        given: dict[str, Any] = {}
        if params.get("expression"):
            given["expression"] = params["expression"]
        if params.get("numerator") or params.get("num"):
            given["numerator"] = params.get("numerator", params.get("num"))
        if params.get("denominator") or params.get("den"):
            given["denominator"] = params.get("denominator", params.get("den"))
        return given

    def get_input_template(self) -> dict:
        return {
            "expression": {
                "type": "str",
                "required": False,
                "description": "F(s) as sympy expression string (e.g. '1/(s**2 + 4*s + 9)')",
            },
            "numerator": {
                "type": "list[float]",
                "required": False,
                "description": "Numerator coefficients [high ... low] (e.g. [1] for constant 1)",
            },
            "denominator": {
                "type": "list[float]",
                "required": False,
                "description": "Denominator coefficients [high ... low] (e.g. [1, 4, 9] for s²+4s+9)",
            },
        }
