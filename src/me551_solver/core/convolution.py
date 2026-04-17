"""Convolution integral (Duhamel's integral) 솔버.

임의 가진 f(t)에 대한 2차/1차 시스템의 시간 응답을 convolution integral로 구한다.

  x(t) = integral_0^t  h(t - tau) f(tau) dtau

여기서 h(t)는 시스템의 impulse response function.

지원 시스템:
  - 2차: mx'' + cx' + kx = f(t)  (underdamped, overdamped, critically damped)
  - 1차: cx' + kx = f(t)

지원 가진:
  - piecewise: 구간별 함수 정의 (예: F0*e^{-at} for 0<t<t1, 0 otherwise)
  - step, ramp, impulse, exponential, harmonic, polynomial 등
"""

from __future__ import annotations

from typing import Any

import numpy as np
import sympy as sp
from sympy import (
    Symbol, Piecewise, exp, sin, cos, sqrt, pi, oo,
    integrate, simplify, nsimplify, Heaviside, DiracDelta,
    Function, Rational, latex, atan2, cancel,
    expand_trig, collect,
)

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Symbolic variables
# ---------------------------------------------------------------------------

t = Symbol("t", real=True, positive=True)
tau = Symbol("tau", real=True, positive=True)


# ---------------------------------------------------------------------------
# Impulse response functions
# ---------------------------------------------------------------------------

def _impulse_response_2nd_order(
    m: sp.Expr, c: sp.Expr, k: sp.Expr, t_var: Symbol = t,
) -> tuple[sp.Expr, dict[str, sp.Expr]]:
    """2차 시스템 mx'' + cx' + kx = f(t)의 impulse response h(t).

    Returns (h(t), info_dict).
    """
    omega_n = sqrt(k / m)
    zeta = c / (2 * sqrt(m * k))

    info: dict[str, sp.Expr] = {
        "omega_n": omega_n,
        "zeta": zeta,
    }

    # Underdamped: zeta < 1
    omega_d = omega_n * sqrt(1 - zeta**2)
    h_under = (1 / (m * omega_d)) * exp(-zeta * omega_n * t_var) * sin(omega_d * t_var)

    # Overdamped: zeta > 1
    s1 = -zeta * omega_n + omega_n * sqrt(zeta**2 - 1)
    s2 = -zeta * omega_n - omega_n * sqrt(zeta**2 - 1)
    h_over = (1 / (m * (s1 - s2))) * (exp(s1 * t_var) - exp(s2 * t_var))

    # Critically damped: zeta = 1
    h_crit = (1 / m) * t_var * exp(-omega_n * t_var)

    info["omega_d"] = omega_d
    info["h_underdamped"] = simplify(h_under)
    info["h_overdamped"] = simplify(h_over)
    info["h_critical"] = simplify(h_crit)

    # Check if zeta is numeric
    zeta_val = None
    try:
        zeta_val = float(zeta)
    except (TypeError, ValueError):
        pass

    if zeta_val is not None:
        if abs(zeta_val - 1.0) < 1e-12:
            h = h_crit
            info["damping_type"] = "critically_damped"
        elif zeta_val < 1.0:
            h = h_under
            info["damping_type"] = "underdamped"
        else:
            h = h_over
            info["damping_type"] = "overdamped"
    else:
        # symbolic zeta -> assume underdamped (most common exam case)
        h = h_under
        info["damping_type"] = "underdamped (assumed)"

    info["h_t"] = simplify(h)
    return h, info


def _impulse_response_1st_order(
    c_val: sp.Expr, k_val: sp.Expr, t_var: Symbol = t,
) -> tuple[sp.Expr, dict[str, sp.Expr]]:
    """1차 시스템 cx' + kx = f(t)의 impulse response h(t).

    h(t) = (1/c) * exp(-(k/c)*t)
    """
    tau_const = c_val / k_val
    h = (1 / c_val) * exp(-t_var / tau_const)

    info: dict[str, sp.Expr] = {
        "time_constant": tau_const,
        "h_t": simplify(h),
    }
    return h, info


# ---------------------------------------------------------------------------
# Piecewise function parser
# ---------------------------------------------------------------------------

def parse_piecewise_input(pieces: list[dict]) -> sp.Expr:
    """구간별 함수 정의를 sympy Piecewise로 변환.

    pieces: list of {"expr": str, "start": str/float, "end": str/float}
        예: [{"expr": "F0*exp(-a*t)", "start": 0, "end": "t1"},
             {"expr": "0", "start": "t1", "end": "inf"}]

    expr 내에서 사용 가능한 변수: t, tau
    사전 정의 심볼: F0, f0, a, b, t1, t2, T, omega, w
    """
    # Pre-defined symbols for convenience
    local_syms = {
        "t": tau,  # convolution에서 적분변수는 tau
        "F0": Symbol("F_0", positive=True),
        "f0": Symbol("f_0", positive=True),
        "f_0": Symbol("f_0", positive=True),
        "a": Symbol("a", positive=True),
        "b": Symbol("b", real=True),
        "alpha": Symbol("alpha", positive=True),
        "t1": Symbol("t_1", positive=True),
        "t2": Symbol("t_2", positive=True),
        "T": Symbol("T", positive=True),
        "omega": Symbol("omega", positive=True),
        "w": Symbol("omega", positive=True),
        "pi": pi,
        "tau": tau,
        "exp": exp,
        "sin": sin,
        "cos": cos,
        "sqrt": sqrt,
    }

    pw_args = []
    for piece in pieces:
        expr_str = str(piece["expr"]).strip()
        start = piece.get("start", 0)
        end = piece.get("end", "inf")

        # Parse expression
        expr = sp.sympify(expr_str, locals=local_syms)

        # Parse bounds
        if isinstance(start, str):
            start_val = sp.sympify(start, locals=local_syms)
        else:
            start_val = nsimplify(start)

        if isinstance(end, str) and end.lower() == "inf":
            end_val = oo
        elif isinstance(end, str):
            end_val = sp.sympify(end, locals=local_syms)
        else:
            end_val = nsimplify(end)

        # Build condition
        if end_val == oo:
            cond = tau >= start_val
        else:
            cond = (tau >= start_val) & (tau < end_val)

        pw_args.append((expr, cond))

    # Default: 0
    pw_args.append((sp.S.Zero, True))
    return Piecewise(*pw_args)


def _build_piecewise_from_params(params: dict) -> tuple[sp.Expr, list[dict], dict]:
    """params에서 piecewise 함수를 구성한다.

    Returns (f_tau_piecewise, pieces_info, symbols_used).
    """
    pieces = params.get("pieces", [])
    if not pieces:
        raise ValueError("'pieces' 리스트가 필요합니다.")

    f_tau = parse_piecewise_input(pieces)
    return f_tau, pieces, {}


# ---------------------------------------------------------------------------
# Preset forcing functions
# ---------------------------------------------------------------------------

def _preset_forcing(params: dict) -> tuple[sp.Expr, str, list[dict]]:
    """프리셋 가진 함수를 생성.

    Returns (f_tau, description, pieces_for_display).
    """
    force_type = params.get("force_type", "step")

    F0 = Symbol("F_0", positive=True) if isinstance(params.get("F0"), str) else nsimplify(params.get("F0", 1))
    a_sym = Symbol("a", positive=True) if isinstance(params.get("a"), str) else nsimplify(params.get("a", 1))
    t1 = Symbol("t_1", positive=True) if isinstance(params.get("t1"), str) else nsimplify(params.get("t1", 1))
    omega = Symbol("omega", positive=True) if isinstance(params.get("omega"), str) else nsimplify(params.get("omega", 1))

    if force_type == "step":
        f = F0
        desc = f"f(t) = {F0} · u(t)"
        pieces = [{"expr": str(F0), "start": 0, "end": "inf"}]

    elif force_type == "ramp":
        f = F0 * tau
        desc = f"f(t) = {F0} · t"
        pieces = [{"expr": f"{F0}*t", "start": 0, "end": "inf"}]

    elif force_type == "exponential":
        f = F0 * exp(-a_sym * tau)
        desc = f"f(t) = {F0} · exp(-{a_sym}·t)"
        pieces = [{"expr": f"{F0}*exp(-{a_sym}*t)", "start": 0, "end": "inf"}]

    elif force_type == "exp_truncated":
        # F0*e^{-at} for 0 < t < t1, 0 for t > t1
        f = Piecewise(
            (F0 * exp(-a_sym * tau), (tau >= 0) & (tau < t1)),
            (sp.S.Zero, True),
        )
        desc = f"f(t) = {{ {F0}·e^(-{a_sym}·t),  0 < t < {t1}\n       {{ 0,              t > {t1}"
        pieces = [
            {"expr": f"{F0}*exp(-{a_sym}*t)", "start": 0, "end": str(t1)},
            {"expr": "0", "start": str(t1), "end": "inf"},
        ]

    elif force_type == "harmonic":
        func = params.get("func", "sin")
        if func == "cos":
            f = F0 * cos(omega * tau)
        else:
            f = F0 * sin(omega * tau)
        desc = f"f(t) = {F0} · {func}({omega}·t)"
        pieces = [{"expr": f"{F0}*{func}({omega}*t)", "start": 0, "end": "inf"}]

    elif force_type == "rectangular_pulse":
        # F0 for 0 < t < t1, 0 for t > t1
        f = Piecewise(
            (F0, (tau >= 0) & (tau < t1)),
            (sp.S.Zero, True),
        )
        desc = f"f(t) = {{ {F0},  0 < t < {t1}\n       {{ 0,   t > {t1}"
        pieces = [
            {"expr": str(F0), "start": 0, "end": str(t1)},
            {"expr": "0", "start": str(t1), "end": "inf"},
        ]

    else:
        raise ValueError(f"Unknown force_type: {force_type}")

    return f, desc, pieces


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class ConvolutionSolver(BaseSolver):
    """Convolution integral (Duhamel's integral) 솔버.

    x(t) = ∫₀ᵗ h(t-τ) f(τ) dτ  (zero initial conditions)
    """

    def solve(self, params: dict) -> SolverResult:
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # -----------------------------------------------------------
        # 0. 시스템 파라미터 파싱
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
        # 1. Impulse response h(t)
        # -----------------------------------------------------------
        if system_order == 2:
            h_t, h_info = _impulse_response_2nd_order(m_val, c_val, k_val, t)
            omega_n = h_info["omega_n"]
            zeta = h_info["zeta"]

            step1_detail = (
                f"시스템: m·x'' + c·x' + k·x = f(t)\n"
                f"ωₙ = √(k/m) = {simplify(omega_n)}\n"
                f"ζ = c/(2√(mk)) = {simplify(zeta)}\n"
                f"Damping type: {h_info['damping_type']}\n"
            )
            if "underdamped" in h_info["damping_type"]:
                omega_d = h_info["omega_d"]
                step1_detail += (
                    f"ωd = ωₙ√(1-ζ²) = {simplify(omega_d)}\n\n"
                    f"h(t) = (1/(m·ωd)) · e^(-ζωₙt) · sin(ωd·t)\n"
                    f"     = {h_info['h_t']}"
                )
            elif h_info["damping_type"] == "critically_damped":
                step1_detail += (
                    f"\nh(t) = (1/m) · t · e^(-ωₙt)\n"
                    f"     = {h_info['h_t']}"
                )
            else:
                step1_detail += f"\nh(t) = {h_info['h_t']}"

            final_answer["omega_n"] = str(simplify(omega_n))
            final_answer["zeta"] = str(simplify(zeta))
            if "underdamped" in h_info["damping_type"]:
                final_answer["omega_d"] = str(simplify(h_info["omega_d"]))
            final_answer["damping_type"] = h_info["damping_type"]
        else:
            h_t, h_info = _impulse_response_1st_order(c_val, k_val, t)
            tau_const = h_info["time_constant"]

            step1_detail = (
                f"시스템: c·x' + k·x = f(t)\n"
                f"시정수 τ = c/k = {simplify(tau_const)}\n\n"
                f"h(t) = (1/c) · e^(-t/τ)\n"
                f"     = {h_info['h_t']}"
            )

        steps.append(("Step 1: Impulse response h(t)", step1_detail))
        final_answer["h_t"] = str(h_info["h_t"])

        # -----------------------------------------------------------
        # 2. 가진 함수 f(t) 정의
        # -----------------------------------------------------------
        input_mode = params.get("input_mode", "preset")

        if input_mode == "piecewise":
            f_tau, pieces, _ = _build_piecewise_from_params(params)
            pieces_desc = []
            for p in params["pieces"]:
                pieces_desc.append(f"  {p['expr']},  {p['start']} < t < {p['end']}")
            f_desc = "f(t) = {\n" + "\n".join(pieces_desc) + "\n}"
        else:
            f_tau, f_desc, pieces = _preset_forcing(params)

        steps.append(("Step 2: 가진 함수 f(t)", f_desc))
        given["f_t"] = f_desc
        final_answer["f_t"] = f_desc

        # -----------------------------------------------------------
        # 3. Convolution integral 설정
        # -----------------------------------------------------------
        # h(t-tau): t를 (t-tau)로 치환
        h_t_minus_tau = h_t.subs(t, t - tau)

        integrand = h_t_minus_tau * f_tau
        integrand_simplified = sp.expand(integrand)

        steps.append((
            "Step 3: Convolution integral 설정",
            f"x(t) = ∫₀ᵗ h(t-τ)·f(τ) dτ\n\n"
            f"h(t-τ) = {simplify(h_t_minus_tau)}\n\n"
            f"피적분함수 = h(t-τ)·f(τ)\n"
            f"         = {integrand_simplified}"
        ))

        # -----------------------------------------------------------
        # 4. 구간별 적분 수행
        # -----------------------------------------------------------
        # Piecewise인 경우 구간별로 분리하여 적분
        if isinstance(f_tau, Piecewise):
            x_t_pieces = self._integrate_piecewise(
                h_t, f_tau, params, steps, final_answer
            )
        else:
            # 단일 구간: 직접 적분
            steps.append(("Step 4: 적분 수행", "단일 구간 적분..."))
            x_t_pieces = self._integrate_single(
                h_t_minus_tau, f_tau, steps, final_answer
            )

        # -----------------------------------------------------------
        # 5. Sanity checks
        # -----------------------------------------------------------
        sanity_parts = []

        # x(0) = 0 확인 (zero IC)
        try:
            x_at_0 = x_t_pieces[0]["x_t"].subs(t, sp.S.Zero)
            x_at_0_simplified = simplify(x_at_0)
            sanity_parts.append(
                f"x(0) = {x_at_0_simplified} "
                f"{'✓' if x_at_0_simplified == 0 else '(should be 0 for zero IC)'}"
            )
        except Exception:
            pass

        # 차원 검증
        if system_order == 2:
            sanity_parts.append(
                f"System: {h_info['damping_type']}, "
                f"ωₙ = {simplify(omega_n)}, ζ = {simplify(zeta)}"
            )
        sanity_parts.append("Zero initial conditions: x(0) = 0, x'(0) = 0")

        return SolverResult(
            problem_type="Convolution Integral (Duhamel's Integral)",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    # -------------------------------------------------------------------
    # Helper: parse symbolic or numeric parameter
    # -------------------------------------------------------------------
    @staticmethod
    def _parse_sym(val, name: str, **assumptions) -> sp.Expr:
        if isinstance(val, str):
            return Symbol(name, **assumptions)
        return nsimplify(val, rational=True)

    # -------------------------------------------------------------------
    # Smart integration: expand trig, factor out t-only parts
    # -------------------------------------------------------------------
    @staticmethod
    def _smart_integrate(integrand: sp.Expr, var: Symbol,
                         lo: sp.Expr, hi: sp.Expr) -> sp.Expr:
        """삼각함수를 전개하고 t-의존 인수를 밖으로 빼서 적분한다.

        핵심 전략:
        1. expand_trig → sin(ωd(t-τ)) = sin(ωd·t)cos(ωd·τ) - cos(ωd·t)sin(ωd·τ)
        2. 전개된 각 항에서 var(τ)를 포함하지 않는 인수를 적분 밖으로 분리
        3. var만 포함하는 부분만 적분 → 표준 ∫e^{ατ}cos(βτ)dτ 형태
        """
        # Step 1: trig 전개 + 곱셈 전개
        expanded = sp.expand(expand_trig(sp.expand(integrand)))

        # Step 2: 합의 각 항을 분리
        terms = sp.Add.make_args(expanded)

        result = sp.S.Zero

        for term in terms:
            # 곱의 인수들을 var-의존/비의존으로 분류
            factors = sp.Mul.make_args(term)
            outer = sp.S.One   # var를 안 포함 (t만 있는 부분)
            inner = sp.S.One   # var를 포함 (tau 부분)

            for f in factors:
                if f.has(var):
                    inner *= f
                else:
                    outer *= f

            # inner만 적분 (표준 적분)
            try:
                int_result = integrate(inner, (var, lo, hi))
                # Piecewise 결과에서 일반 경우(첫 번째 분기)만 추출
                if isinstance(int_result, Piecewise) and len(int_result.args) >= 1:
                    int_result = int_result.args[0][0]
                if int_result.has(sp.Integral):
                    # 한 번 더 시도: simplify 후 재적분
                    int_result = integrate(
                        sp.simplify(inner), (var, lo, hi)
                    )
                    if isinstance(int_result, Piecewise):
                        int_result = int_result.args[0][0]
            except Exception:
                int_result = sp.Integral(inner, (var, lo, hi))

            result += outer * int_result

        # 최종 정리
        try:
            result = simplify(result)
        except Exception:
            pass

        return result

    # -------------------------------------------------------------------
    # Piecewise integration
    # -------------------------------------------------------------------
    def _integrate_piecewise(
        self, h_t_expr: sp.Expr, f_tau_pw: sp.Piecewise,
        params: dict, steps: list, final_answer: dict,
    ) -> list[dict]:
        """Piecewise f(τ)에 대해 구간별 convolution 적분을 수행."""

        h_t_tau = h_t_expr.subs(t, t - tau)

        # f(τ)의 비-zero 구간과 경계 추출
        piece_intervals = self._extract_piece_intervals(f_tau_pw)

        if not piece_intervals:
            steps.append(("Step 4: 적분", "f(t) = 0 → x(t) = 0"))
            final_answer["x_t"] = "0"
            return [{"region": "all t", "x_t": sp.S.Zero}]

        # 경계점 수집 (0 제외, ∞ 제외)
        bps = sorted(set(
            b for _, lo, hi in piece_intervals
            for b in [lo, hi] if b != 0 and b != oo and not b.is_infinite
        ), key=lambda x: x)

        step_num = 4

        if not bps:
            # 경계점 없음 (step 등): 단일 구간
            expr_0 = piece_intervals[0][0]
            integrand = h_t_tau * expr_0
            integral_result = self._smart_integrate(integrand, tau, sp.S.Zero, t)

            steps.append((
                f"Step {step_num}: 적분 수행 (t > 0)",
                f"x(t) = ∫₀ᵗ h(t-τ)·f(τ) dτ\n"
                f"     = ∫₀ᵗ ({simplify(h_t_tau)})·({expr_0}) dτ\n\n"
                f"     = {integral_result}"
            ))
            final_answer["x_t"] = str(integral_result)
            return [{"region": "t > 0", "x_t": integral_result}]

        # 경계점이 있는 경우: 시간 영역별 적분
        time_regions = []

        for region_idx in range(len(bps) + 1):
            if region_idx == 0:
                region_label = f"0 < t < {bps[0]}"
            elif region_idx < len(bps):
                region_label = f"{bps[region_idx-1]} < t < {bps[region_idx]}"
            else:
                region_label = f"t > {bps[-1]}"

            total_integral = sp.S.Zero
            integral_parts = []

            for expr_i, lo_i, hi_i in piece_intervals:
                # 적분 하한: max(lo_i, 0) = lo_i (항상 >= 0)
                int_lo = lo_i

                # 적분 상한: min(hi_i, t)
                if hi_i == oo or hi_i.is_infinite:
                    int_hi = t
                else:
                    if region_idx == 0:
                        # t < bps[0] <= hi_i → 상한 = t
                        int_hi = t
                    else:
                        # t > hi_i → 상한 = hi_i (f가 끝난 구간)
                        int_hi = hi_i

                integrand = h_t_tau * expr_i
                result = self._smart_integrate(integrand, tau, int_lo, int_hi)

                if result != 0:
                    total_integral += result
                    integral_parts.append(
                        f"∫_{{{int_lo}}}^{{{int_hi}}} h(t-τ)·({expr_i}) dτ = {result}"
                    )

            total_integral = simplify(total_integral)

            step_num += 1
            detail = f"Region: {region_label}\n\n"
            if integral_parts:
                detail += "\n".join(integral_parts)
                detail += f"\n\nx(t) = {total_integral}"
            else:
                detail += "x(t) = 0"

            steps.append((f"Step {step_num}: {region_label}", detail))
            time_regions.append({"region": region_label, "x_t": total_integral})

        # Final summary
        final_answer["x_t_regions"] = [
            {"region": r["region"], "x_t": str(r["x_t"])}
            for r in time_regions
        ]
        summary_lines = [f"  {r['region']}:  x(t) = {r['x_t']}" for r in time_regions]
        steps.append((
            f"Step {step_num + 1}: 최종 응답 x(t)",
            "x(t) = {\n" + "\n".join(summary_lines) + "\n}"
        ))

        return time_regions

    # -------------------------------------------------------------------
    # Single-interval integration
    # -------------------------------------------------------------------
    def _integrate_single(
        self, h_t_tau: sp.Expr, f_tau_expr: sp.Expr,
        steps: list, final_answer: dict,
    ) -> list[dict]:
        integrand = h_t_tau * f_tau_expr
        result = self._smart_integrate(integrand, tau, sp.S.Zero, t)

        steps.append((
            "Step 4: 적분 수행",
            f"x(t) = ∫₀ᵗ ({simplify(h_t_tau)})·({f_tau_expr}) dτ\n\n"
            f"     = {result}"
        ))
        final_answer["x_t"] = str(result)
        return [{"region": "t > 0", "x_t": result}]

    # -------------------------------------------------------------------
    # Piecewise interval extraction
    # -------------------------------------------------------------------
    @staticmethod
    def _extract_piece_intervals(pw: sp.Piecewise) -> list[tuple]:
        """Piecewise에서 (expr, lo, hi) 리스트를 추출.

        SymPy가 조건을 다양하게 단순화하므로 여러 패턴을 처리:
          - And(tau >= a, tau < b) → (a, b)
          - StrictGreaterThan(b, tau) i.e. b > tau → (0, b)  (tau positive)
          - StrictLessThan(tau, b) i.e. tau < b → (0, b)
          - tau >= a → (a, ∞)
          - True → default (skip)
        """
        intervals = []
        for expr_i, cond_i in pw.args:
            if expr_i == 0 or cond_i is sp.true:
                continue

            lo, hi = sp.S.Zero, oo

            # Case 1: And(rel1, rel2)
            if isinstance(cond_i, sp.And):
                for rel in cond_i.args:
                    bound = ConvolutionSolver._parse_relational(rel)
                    if bound:
                        if bound[0] == "lo":
                            lo = bound[1]
                        else:
                            hi = bound[1]
            # Case 2: single relational (e.g. t_1 > tau, tau < t_1)
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
        """단일 관계식에서 tau의 상한/하한을 추출.

        Returns ("lo", value) or ("hi", value) or None.
        """
        if not hasattr(rel, 'lhs'):
            return None

        lhs, rhs = rel.lhs, rel.rhs

        # tau >= a  or  tau > a  → lo = a
        if lhs == tau and isinstance(rel, (sp.core.relational.GreaterThan,
                                           sp.core.relational.StrictGreaterThan)):
            return ("lo", rhs)

        # a <= tau  or  a < tau  → lo = a
        if rhs == tau and isinstance(rel, (sp.core.relational.LessThan,
                                           sp.core.relational.StrictLessThan)):
            return ("lo", lhs)

        # tau <= b  or  tau < b  → hi = b
        if lhs == tau and isinstance(rel, (sp.core.relational.LessThan,
                                           sp.core.relational.StrictLessThan)):
            return ("hi", rhs)

        # b >= tau  or  b > tau  → hi = b
        if rhs == tau and isinstance(rel, (sp.core.relational.GreaterThan,
                                           sp.core.relational.StrictGreaterThan)):
            return ("hi", lhs)

        return None

    def get_input_template(self) -> dict:
        return {
            "system_order": {
                "type": "int",
                "default": 2,
                "description": "시스템 차수 (1 or 2)",
            },
            "m": {"type": "float/str", "description": "질량 (2차 시스템)"},
            "c": {"type": "float/str", "description": "감쇠 계수"},
            "k": {"type": "float/str", "description": "강성 계수"},
            "input_mode": {
                "type": "str",
                "description": "'preset' (프리셋) 또는 'piecewise' (구간별 정의)",
            },
            "force_type": {
                "type": "str",
                "description": "프리셋: 'step', 'ramp', 'exponential', "
                               "'exp_truncated', 'harmonic', 'rectangular_pulse'",
            },
            "pieces": {
                "type": "list[dict]",
                "description": "piecewise 모드: [{'expr': 'F0*exp(-a*t)', 'start': 0, 'end': 't1'}, ...]",
            },
        }
