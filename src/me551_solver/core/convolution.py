"""Convolution integral (Duhamel's integral) 솔버.

임의 가진 f(t)에 대한 2차/1차 시스템의 시간 응답을 convolution integral로 구한다.

  x(t) = integral_0^t  h(t - tau) f(tau) dtau

여기서 h(t)는 시스템의 impulse response function.

핵심 적분 공식 (SymPy integrate 미사용, 직접 해석해):
  ∫₀ᴸ e^{σ̃τ} sin(β(t-τ)) dτ
  = [e^{σ̃L}(σ̃·sin(β(t-L)) + β·cos(β(t-L))) - σ̃·sin(βt) - β·cos(βt)] / (σ̃²+β²)
"""

from __future__ import annotations

from typing import Any

import numpy as np
import sympy as sp
from sympy import (
    Symbol, Piecewise, exp, sin, cos, sqrt, pi, oo,
    integrate, simplify, nsimplify, Heaviside,
    Rational, atan2, cancel,
)

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def _pretty(expr: sp.Expr) -> str:
    """SymPy 식을 교과서 표기 문자열로 변환."""
    s = str(expr)
    return (s
            .replace("zeta_omega_n", "ζωₙ")
            .replace("omega_d", "ωd")
            .replace("F_0", "F₀")
            .replace("t_1", "t₁")
            .replace("t_2", "t₂")
            )


# ---------------------------------------------------------------------------
# Symbolic variables
# ---------------------------------------------------------------------------

t = Symbol("t", real=True, positive=True)
tau = Symbol("tau", real=True, positive=True)

# 표준 표기 심볼
_zwn = Symbol("zeta_omega_n", positive=True)   # ζωₙ
_wd = Symbol("omega_d", positive=True)          # ωd
_m_s = Symbol("m", positive=True)


# ---------------------------------------------------------------------------
# Analytical convolution formula (NO SymPy integrate)
# ---------------------------------------------------------------------------

def _conv_exp_sin(sigma: sp.Expr, beta: sp.Expr, L: sp.Expr,
                  t_var: Symbol = t) -> sp.Expr:
    """∫₀ᴸ e^{σ·τ} sin(β(t-τ)) dτ 의 해석해.

    표준 공식:
      = [e^{σL}(σ·sin(β(t-L)) + β·cos(β(t-L))) - σ·sin(βt) - β·cos(βt)]
        / (σ² + β²)

    증명: sin(β(t-τ)) = sin(βt)cos(βτ) - cos(βt)sin(βτ) 전개 후
         ∫e^{στ}cos(βτ)dτ, ∫e^{στ}sin(βτ)dτ 표준 적분 적용.
    """
    D = sigma**2 + beta**2  # 분모
    return (
        exp(sigma * L) * (sigma * sin(beta * (t_var - L)) + beta * cos(beta * (t_var - L)))
        - sigma * sin(beta * t_var) - beta * cos(beta * t_var)
    ) / D


# ---------------------------------------------------------------------------
# Impulse response
# ---------------------------------------------------------------------------

def _impulse_response_2nd(m, c, k):
    """2차 시스템 파라미터 추출. Returns info dict."""
    omega_n = sqrt(k / m)
    zeta = c / (2 * sqrt(m * k))
    omega_d = omega_n * sqrt(1 - zeta**2)
    zeta_omega_n = zeta * omega_n  # = c/(2m)

    zeta_val = None
    try:
        zeta_val = float(zeta)
    except (TypeError, ValueError):
        pass

    if zeta_val is not None:
        if abs(zeta_val - 1.0) < 1e-12:
            dtype = "critically_damped"
        elif zeta_val > 1.0:
            dtype = "overdamped"
        else:
            dtype = "underdamped"
    else:
        dtype = "underdamped (assumed)"

    return {
        "omega_n": omega_n,
        "zeta": zeta,
        "omega_d_expr": simplify(omega_d),
        "zeta_omega_n_expr": simplify(zeta_omega_n),
        "damping_type": dtype,
    }


# ---------------------------------------------------------------------------
# Preset forcing: extract (amplitude, exp_rate) for each piece
# ---------------------------------------------------------------------------

def _parse_forcing(params: dict):
    """가진 함수를 구간별 (A, p, lo, hi) 리스트로 파싱.

    f(τ) = A·e^{p·τ} on [lo, hi)
    step:  A=F0, p=0
    exp:   A=F0, p=-a
    """
    force_type = params.get("force_type", "step")

    F0 = Symbol("F_0", positive=True) if isinstance(params.get("F0"), str) else nsimplify(params.get("F0", 1))
    a = Symbol("a", positive=True) if isinstance(params.get("a"), str) else nsimplify(params.get("a", 1))
    t1 = Symbol("t_1", positive=True) if isinstance(params.get("t1"), str) else nsimplify(params.get("t1", 1))

    if force_type == "step":
        desc = f"f(t) = {F0} · u(t)"
        pieces = [(F0, sp.S.Zero, sp.S.Zero, oo)]  # (A, p, lo, hi)

    elif force_type == "exponential":
        desc = f"f(t) = {F0} · exp(-{a}·t)"
        pieces = [(F0, -a, sp.S.Zero, oo)]

    elif force_type == "exp_truncated":
        desc = (
            f"f(t) = {{ {F0}·exp(-{a}·t),  0 < t < {t1}\n"
            f"       {{ 0,              t > {t1}"
        )
        pieces = [(F0, -a, sp.S.Zero, t1)]

    elif force_type == "rectangular_pulse":
        desc = (
            f"f(t) = {{ {F0},  0 < t < {t1}\n"
            f"       {{ 0,   t > {t1}"
        )
        pieces = [(F0, sp.S.Zero, sp.S.Zero, t1)]

    elif force_type == "ramp":
        desc = f"f(t) = {F0} · t"
        # ramp은 exp 형태가 아니므로 특수 처리 플래그
        pieces = [("ramp", F0, sp.S.Zero, oo)]

    elif force_type == "harmonic":
        func = params.get("func", "sin")
        omega = Symbol("omega", positive=True) if isinstance(params.get("omega"), str) else nsimplify(params.get("omega", 1))
        desc = f"f(t) = {F0} · {func}({omega}·t)"
        pieces = [("harmonic", F0, func, omega, sp.S.Zero, oo)]

    else:
        raise ValueError(f"Unknown force_type: {force_type}")

    return desc, pieces


def _normalize_expr_str(s: str) -> str:
    """사용자 입력을 SymPy sympify가 파싱 가능한 형태로 정규화.

    e^(...)  → exp(...)
    e**(...) → exp(...)
    E^(...)  → exp(...)
    ^        → **
    """
    import re
    # e^(...) 또는 e**(...) → exp(...)
    s = re.sub(r'\be\s*\^\s*\(', 'exp(', s)
    s = re.sub(r'\be\s*\*\*\s*\(', 'exp(', s)
    s = re.sub(r'\bE\s*\^\s*\(', 'exp(', s)
    # 남은 ^ → **
    s = s.replace('^', '**')
    return s


def _parse_piecewise_forcing(params: dict):
    """사용자 정의 구간별 입력을 (A, p, lo, hi) 형태로 변환 시도."""
    raw_pieces = params.get("pieces", [])
    desc_lines = []
    parsed = []

    local_syms = {
        "t": Symbol("_tau_dummy", positive=True),
        "F0": Symbol("F_0", positive=True),
        "f0": Symbol("f_0", positive=True),
        "f_0": Symbol("f_0", positive=True),
        "a": Symbol("a", positive=True),
        "b": Symbol("b", real=True),
        "t1": Symbol("t_1", positive=True),
        "t2": Symbol("t_2", positive=True),
        "T": Symbol("T", positive=True),
        "omega": Symbol("omega", positive=True),
        "w": Symbol("omega", positive=True),
        "exp": exp, "sin": sin, "cos": cos, "sqrt": sqrt, "pi": pi,
    }
    _tau_d = local_syms["t"]

    for p in raw_pieces:
        expr_str = _normalize_expr_str(str(p["expr"]).strip())
        start = p.get("start", "0")
        end = p.get("end", "inf")

        lo_str = _normalize_expr_str(str(start))
        hi_str = _normalize_expr_str(str(end))

        lo = sp.sympify(lo_str, locals=local_syms) if lo_str != "0" else sp.S.Zero
        hi = oo if hi_str.lower() == "inf" else sp.sympify(hi_str, locals=local_syms)

        desc_lines.append(f"  {p['expr']},  {start} < t < {end}")

        expr = sp.sympify(expr_str, locals=local_syms)

        if expr == 0:
            continue

        # A·exp(p·τ) 패턴 감지
        if not expr.has(_tau_d):
            # 상수 → f(τ) = A, p = 0
            parsed.append((expr, sp.S.Zero, lo, hi))
        else:
            # A*exp(p*t) 패턴 매칭
            A_w = sp.Wild("A", exclude=[_tau_d])
            p_w = sp.Wild("p", exclude=[_tau_d])
            m_exp = expr.match(A_w * exp(p_w * _tau_d))
            if m_exp:
                parsed.append((m_exp[A_w], m_exp[p_w], lo, hi))
            else:
                # exp 없는 곱? A*t 같은 경우 체크
                # 마지막으로 A*exp(stuff) 형태 시도 (stuff가 단순 곱이 아닌 경우)
                # 예: F0*exp(-a*t) 에서 -a*t → p = -a
                atoms_exp = list(expr.atoms(sp.exp))
                if len(atoms_exp) == 1:
                    # exp 하나 포함 → 계수와 지수 분리
                    exp_part = atoms_exp[0]
                    coeff = expr / exp_part
                    coeff = simplify(coeff)
                    # 지수에서 t 계수 추출
                    exponent = exp_part.args[0]
                    p_val = exponent.coeff(_tau_d)
                    if p_val != 0 and simplify(exponent - p_val * _tau_d) == 0:
                        parsed.append((coeff, p_val, lo, hi))
                    else:
                        # fallback
                        parsed.append(("general", expr.subs(_tau_d, tau), lo, hi))
                else:
                    parsed.append(("general", expr.subs(_tau_d, tau), lo, hi))

    desc = "f(t) = {\n" + "\n".join(desc_lines) + "\n}"
    return desc, parsed


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class ConvolutionSolver(BaseSolver):
    """Convolution integral (Duhamel's integral) 솔버."""

    def solve(self, params: dict) -> SolverResult:
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # -----------------------------------------------------------
        # 0. 시스템 파라미터
        # -----------------------------------------------------------
        system_order = params.get("system_order", 2)

        m_val = self._parse_sym(params.get("m", 1), "m", positive=True) if system_order >= 2 else None
        c_val = self._parse_sym(params.get("c", 1), "c", positive=True)
        k_val = self._parse_sym(params.get("k", 1), "k", positive=True)

        given: dict[str, Any] = {"system_order": system_order}
        if system_order == 2:
            given["equation"] = f"m·x'' + c·x' + k·x = f(t),  m={m_val}, c={c_val}, k={k_val}"
        else:
            given["equation"] = f"c·x' + k·x = f(t),  c={c_val}, k={k_val}"

        # -----------------------------------------------------------
        # 1. Impulse response
        # -----------------------------------------------------------
        use_analytic = (system_order == 2)

        if use_analytic:
            info = _impulse_response_2nd(m_val, c_val, k_val)
            omega_n = info["omega_n"]
            zeta = info["zeta"]
            dtype = info["damping_type"]

            step1 = (
                f"시스템: m·x'' + c·x' + k·x = f(t)\n"
                f"ωₙ = √(k/m) = {_pretty(simplify(omega_n))}\n"
                f"ζ = c/(2√(mk)) = {_pretty(simplify(zeta))}\n"
                f"Damping type: {dtype}\n"
                f"ωd = ωₙ√(1-ζ²) = {_pretty(info['omega_d_expr'])}\n"
                f"ζωₙ = c/(2m) = {_pretty(info['zeta_omega_n_expr'])}\n\n"
                f"h(t) = (1/(m·ωd)) · exp(-ζωₙ·t) · sin(ωd·t)"
            )
            final_answer["omega_n"] = _pretty(simplify(omega_n))
            final_answer["zeta"] = _pretty(simplify(zeta))
            final_answer["omega_d"] = _pretty(info["omega_d_expr"])
            final_answer["damping_type"] = dtype
            final_answer["h_t"] = "exp(-ζωₙ·t)·sin(ωd·t)/(m·ωd)"
        else:
            tau_c = c_val / k_val
            step1 = (
                f"시스템: c·x' + k·x = f(t)\n"
                f"시정수 τ = c/k = {simplify(tau_c)}\n\n"
                f"h(t) = (1/c) · exp(-t/τ)"
            )
            final_answer["h_t"] = f"(1/{c_val})*exp(-{k_val}/{c_val}*t)"

        steps.append(("Step 1: Impulse response h(t)", step1))

        # -----------------------------------------------------------
        # 2. 가진 함수
        # -----------------------------------------------------------
        input_mode = params.get("input_mode", "preset")
        if input_mode == "piecewise":
            f_desc, force_pieces = _parse_piecewise_forcing(params)
        else:
            f_desc, force_pieces = _parse_forcing(params)

        steps.append(("Step 2: 가진 함수 f(t)", _pretty(f_desc)))
        given["f_t"] = _pretty(f_desc)
        final_answer["f_t"] = _pretty(f_desc)

        # -----------------------------------------------------------
        # 3. Convolution integral 설정
        # -----------------------------------------------------------
        steps.append((
            "Step 3: Convolution integral 설정",
            "x(t) = (1/(m·ωd)) ∫₀ᵗ exp(-ζωₙ(t-τ))·sin(ωd(t-τ))·f(τ) dτ\n\n"
            "f(τ) = A·exp(p·τ) 형태일 때:\n"
            "  피적분함수 = (A/(m·ωd))·exp(-ζωₙ·t)·exp((ζωₙ+p)·τ)·sin(ωd(t-τ))\n\n"
            "σ̃ = ζωₙ + p 로 치환하면 표준 적분:\n"
            "  ∫₀ᴸ exp(σ̃τ)·sin(ωd(t-τ)) dτ\n"
            "  = [exp(σ̃L)·(σ̃·sin(ωd(t-L)) + ωd·cos(ωd(t-L))) - σ̃·sin(ωd·t) - ωd·cos(ωd·t)]\n"
            "    / (σ̃² + ωd²)"
        ))

        # -----------------------------------------------------------
        # 4. 구간별 적분 (해석적 공식 직접 적용)
        # -----------------------------------------------------------
        # 경계점 수집
        bps = set()
        normal_pieces = []  # (A, p, lo, hi)
        for piece in force_pieces:
            if isinstance(piece[0], str) and piece[0] in ("ramp", "harmonic", "general"):
                # fallback to SymPy integrate
                normal_pieces.append(piece)
                continue
            A, p, lo, hi = piece
            normal_pieces.append((A, p, lo, hi))
            if hi != oo and not hi.is_infinite:
                bps.add(hi)

        bps = sorted(bps)

        time_regions = []
        step_num = 3

        if not bps:
            # 단일 구간 (step, exponential 등)
            piece = normal_pieces[0]
            if isinstance(piece[0], str):
                # fallback
                x_t = self._fallback_integrate(piece, system_order, m_val, c_val, k_val)
            else:
                A, p, lo, hi = piece
                x_t = self._analytic_response(A, p, t, m_val)
            x_t = simplify(x_t)

            step_num += 1
            steps.append((
                f"Step {step_num}: 적분 수행 (t > 0)",
                f"σ̃ = ζωₙ + ({_pretty(p)}) = {_pretty(_zwn + p)}\n"
                f"L = t\n\n"
                f"x(t) = {_pretty(x_t)}" if not isinstance(piece[0], str) else f"x(t) = {_pretty(x_t)}"
            ))
            time_regions.append({"region": "t > 0", "x_t": x_t})
        else:
            # 구간별 적분
            for region_idx in range(len(bps) + 1):
                if region_idx == 0:
                    region_label = f"0 < t < {_pretty(bps[0])}"
                else:
                    region_label = f"t > {_pretty(bps[-1])}"

                total = sp.S.Zero
                detail_lines = []

                for piece in normal_pieces:
                    if isinstance(piece[0], str):
                        continue
                    A, p, lo, hi = piece

                    sigma_tilde = _zwn + p

                    if hi == oo or hi.is_infinite:
                        L = t
                    else:
                        if region_idx == 0:
                            L = t  # t < breakpoint, 적분 상한 = t
                        else:
                            L = hi  # t > breakpoint, 적분 상한 = hi

                    # 해석적 공식 적용
                    conv_integral = _conv_exp_sin(sigma_tilde, _wd, L)

                    # x(t) = A/(m·ωd) · exp(-ζωₙ·t) · conv_integral
                    x_piece = (A / (_m_s * _wd)) * exp(-_zwn * t) * conv_integral
                    x_piece = sp.expand(x_piece)

                    detail_lines.append(
                        f"σ̃ = ζωₙ + ({_pretty(p)}) = {_pretty(sigma_tilde)}\n"
                        f"L = {_pretty(L)}\n"
                        f"∫₀^{{{_pretty(L)}}} exp(σ̃τ)·sin(ωd(t-τ)) dτ\n"
                        f"  = [exp(σ̃·{_pretty(L)})·(σ̃·sin(ωd(t-{_pretty(L)})) + ωd·cos(ωd(t-{_pretty(L)})))\n"
                        f"     - σ̃·sin(ωd·t) - ωd·cos(ωd·t)] / (σ̃² + ωd²)"
                    )

                    total += x_piece

                total = simplify(total)

                step_num += 1
                detail = f"Region: {region_label}\n\n" + "\n\n".join(detail_lines)
                detail += f"\n\nx(t) = {_pretty(total)}"

                steps.append((f"Step {step_num}: {region_label}", detail))
                time_regions.append({"region": region_label, "x_t": total})

        # -----------------------------------------------------------
        # 5. 최종 응답 정리
        # -----------------------------------------------------------
        final_answer["x_t_regions"] = [
            {"region": r["region"], "x_t": _pretty(r["x_t"])}
            for r in time_regions
        ]
        summary_lines = [f"  {r['region']}:  x(t) = {_pretty(r['x_t'])}" for r in time_regions]
        steps.append((
            f"Step {step_num + 1}: 최종 응답 x(t)",
            "x(t) = {\n" + "\n".join(summary_lines) + "\n}"
        ))

        # -----------------------------------------------------------
        # 6. Sanity check
        # -----------------------------------------------------------
        sanity_parts = []
        try:
            x0 = time_regions[0]["x_t"].subs(t, sp.S.Zero)
            x0 = simplify(x0)
            sanity_parts.append(f"x(0) = {x0} {'✓' if x0 == 0 else ''}")
        except Exception:
            pass

        if use_analytic:
            sanity_parts.append(f"System: {dtype}, ωₙ = {_pretty(simplify(omega_n))}, ζ = {_pretty(simplify(zeta))}")
        sanity_parts.append("Zero initial conditions: x(0) = 0, x'(0) = 0")

        return SolverResult(
            problem_type="Convolution Integral (Duhamel's Integral)",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    # -------------------------------------------------------------------
    # Helpers
    # -------------------------------------------------------------------
    @staticmethod
    def _parse_sym(val, name: str, **assumptions) -> sp.Expr:
        if isinstance(val, str):
            return Symbol(name, **assumptions)
        return nsimplify(val, rational=True)

    @staticmethod
    def _analytic_response(A, p, L, m_val):
        """단일 구간 A·e^{pt} 가진에 대한 해석적 응답 (∫₀ᴸ).

        x(t) = A/(m·ωd) · e^{-ζωₙt} · _conv_exp_sin(ζωₙ+p, ωd, L)
        """
        sigma = _zwn + p
        conv = _conv_exp_sin(sigma, _wd, L)
        return (A / (_m_s * _wd)) * exp(-_zwn * t) * conv

    @staticmethod
    def _fallback_integrate(piece, order, m, c, k):
        """비-exponential 가진에 대한 SymPy integrate fallback."""
        # TODO: ramp, harmonic 등
        return sp.S.Zero

    # -------------------------------------------------------------------
    # Piecewise condition parsing (for custom piecewise input)
    # -------------------------------------------------------------------
    @staticmethod
    def _extract_piece_intervals(pw: sp.Piecewise) -> list[tuple]:
        intervals = []
        for expr_i, cond_i in pw.args:
            if expr_i == 0 or cond_i is sp.true:
                continue
            lo, hi = sp.S.Zero, oo
            if isinstance(cond_i, sp.And):
                for rel in cond_i.args:
                    bound = ConvolutionSolver._parse_relational(rel)
                    if bound:
                        if bound[0] == "lo":
                            lo = bound[1]
                        else:
                            hi = bound[1]
            elif isinstance(cond_i, sp.core.relational.Relational):
                bound = ConvolutionSolver._parse_relational(cond_i)
                if bound:
                    if bound[0] == "lo":
                        lo = bound[1]
                    else:
                        hi = bound[1]
            intervals.append((expr_i, lo, hi))
        return intervals

    @staticmethod
    def _parse_relational(rel) -> tuple | None:
        if not hasattr(rel, 'lhs'):
            return None
        lhs, rhs = rel.lhs, rel.rhs
        if lhs == tau and isinstance(rel, (sp.core.relational.GreaterThan, sp.core.relational.StrictGreaterThan)):
            return ("lo", rhs)
        if rhs == tau and isinstance(rel, (sp.core.relational.LessThan, sp.core.relational.StrictLessThan)):
            return ("lo", lhs)
        if lhs == tau and isinstance(rel, (sp.core.relational.LessThan, sp.core.relational.StrictLessThan)):
            return ("hi", rhs)
        if rhs == tau and isinstance(rel, (sp.core.relational.GreaterThan, sp.core.relational.StrictGreaterThan)):
            return ("hi", lhs)
        return None

    def get_input_template(self) -> dict:
        return {
            "system_order": {"type": "int", "default": 2, "description": "시스템 차수 (1 or 2)"},
            "m": {"type": "float/str", "description": "질량 (2차 시스템)"},
            "c": {"type": "float/str", "description": "감쇠 계수"},
            "k": {"type": "float/str", "description": "강성 계수"},
            "input_mode": {"type": "str", "description": "'preset' or 'piecewise'"},
            "force_type": {"type": "str", "description": "'step','ramp','exponential','exp_truncated','harmonic','rectangular_pulse'"},
            "pieces": {"type": "list[dict]", "description": "[{'expr':'F0*exp(-a*t)','start':0,'end':'t1'},...]"},
        }
