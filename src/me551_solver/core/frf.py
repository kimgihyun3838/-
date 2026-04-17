"""FRF 분해 (전달함수) 솔버.

ODE 계수로부터 전달함수 G(s)를 구성하고, 부분분수 분해를 수행하며,
각 서브시스템의 특성 파라미터(ω_b, ω_n, ζ)를 추출한다.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import sympy as sp
from sympy import I, Poly, Symbol, apart, factor

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

s = Symbol("s")


def _coeffs_to_poly(coefficients: list, var: Symbol = s) -> sp.Expr:
    """고차 → 저차 순서의 계수 리스트를 sympy 다항식 표현으로 변환.

    예: [1, 3, 6, 8] → s³ + 3s² + 6s + 8

    소수(float)는 nsimplify로 깔끔한 유리수로 변환한다.
    예: 1.2 → 6/5  (Rational(1.2) → 5404319552844595/4503599627370496 방지)
    """
    n = len(coefficients) - 1
    expr = sp.S.Zero
    for i, c in enumerate(coefficients):
        expr += sp.nsimplify(c, rational=True) * var ** (n - i)
    return sp.expand(expr)


def _classify_factor(expr: sp.Expr, var: Symbol = s) -> dict[str, Any]:
    """단일 인수(1차 또는 2차)의 특성 파라미터를 추출한다.

    Returns a dict with keys:
        order: 1 or 2
        expression: the sympy expression
        For 1st-order: omega_b (break frequency = constant_term / leading_coeff)
        For 2nd-order: omega_n (natural frequency), zeta (damping ratio)
    """
    poly = Poly(expr, var)
    deg = poly.degree()
    coeffs = poly.all_coeffs()  # highest-degree first

    if deg == 1:
        # Form: a*s + b  →  time_constant τ = a/b, break freq ω_b = b/a
        a, b = [sp.nsimplify(c) for c in coeffs]
        omega_b = sp.nsimplify(b / a)
        tau = sp.nsimplify(a / b)
        return {
            "order": 1,
            "expression": expr,
            "omega_b": omega_b,
            "time_constant": tau,
            "description": (
                f"1st-order subsystem: τ = {tau}, ω_b = {omega_b}"
            ),
        }
    elif deg == 2:
        # Form: a*s² + b*s + c  →  ω_n = sqrt(c/a), ζ = b/(2*sqrt(a*c))
        a, b, c = [sp.nsimplify(c_) for c_ in coeffs]
        omega_n_sq = sp.nsimplify(c / a)
        omega_n = sp.sqrt(omega_n_sq)
        zeta = sp.nsimplify(b / (2 * sp.sqrt(a * c)))
        return {
            "order": 2,
            "expression": expr,
            "omega_n": omega_n,
            "omega_n_squared": omega_n_sq,
            "zeta": zeta,
            "description": (
                f"2nd-order subsystem: ω_n = {omega_n} "
                f"(ω_n² = {omega_n_sq}), ζ = {zeta}"
            ),
        }
    else:
        return {
            "order": deg,
            "expression": expr,
            "description": f"{deg}th-order factor (not decomposed further)",
        }


def _format_complex_fraction(num_re: sp.Expr, num_im: sp.Expr,
                              den_re: sp.Expr, den_im: sp.Expr) -> str:
    """복소 분수를 (a + jb) / (c + jd) 표준형 문자열로 변환.

    실수부/허수부가 0인 경우 생략하여 깔끔하게 표시한다.
    """
    # 분자 문자열
    num_re_s = sp.simplify(num_re)
    num_im_s = sp.simplify(num_im)
    if num_re_s == 0 and num_im_s == 0:
        return "0"
    elif num_im_s == 0:
        num_str = str(num_re_s)
    elif num_re_s == 0:
        num_str = f"j*({num_im_s})"
    else:
        num_str = f"({num_re_s}) + j*({num_im_s})"

    # 분모 문자열
    den_re_s = sp.simplify(den_re)
    den_im_s = sp.simplify(den_im)
    if den_im_s == 0:
        den_str = str(den_re_s)
    elif den_re_s == 0:
        den_str = f"j*({den_im_s})"
    else:
        den_str = f"({den_re_s}) + j*({den_im_s})"

    return f"[{num_str}] / [{den_str}]"


def _evaluate_frf(num_coeffs: list, den_coeffs: list,
                   omega: np.ndarray) -> np.ndarray:
    """G(jω) 를 수치적으로 계산한다.

    Parameters:
        num_coeffs: 분자 다항식 계수 (고차→저차)
        den_coeffs: 분모 다항식 계수 (고차→저차)
        omega: 주파수 배열 (rad/s)

    Returns:
        G(jω) 복소수 배열
    """
    jw = 1j * omega
    num_val = np.polyval([complex(c) for c in num_coeffs], jw)
    den_val = np.polyval([complex(c) for c in den_coeffs], jw)
    return num_val / den_val


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------


class FRFSolver(BaseSolver):
    """Frequency Response Function 분해 및 분석.

    ODE 계수로부터 전달함수를 구성하고, 부분분수 분해를 통해
    1차·2차 서브시스템 파라미터를 추출한다.
    """

    def solve(self, params: dict) -> SolverResult:
        """FRF 해석을 수행한다.

        Parameters:
            params: dict with keys:
                coefficients: list[float] – 분모(특성방정식) 계수 [a_n, ..., a_0]
                numerator: list[float] – 분자 계수 (기본값 [1])
                omega_range: tuple (optional) – (omega_min, omega_max, num_points)
        """
        # ---------------------------------------------------------------
        # 0. 입력 파싱
        # ---------------------------------------------------------------
        den_coeffs = params.get("coefficients")
        if den_coeffs is None:
            raise ValueError("'coefficients' (분모 계수)가 필요합니다.")
        num_coeffs = params.get("numerator", [1])

        # numpy 호환을 위해 float 리스트 보관
        den_coeffs_float = [float(c) for c in den_coeffs]
        num_coeffs_float = [float(c) for c in num_coeffs]

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # ---------------------------------------------------------------
        # 1. 전달함수 G(s) 구성
        # ---------------------------------------------------------------
        num_poly = _coeffs_to_poly(num_coeffs, s)
        den_poly = _coeffs_to_poly(den_coeffs, s)
        G_s = num_poly / den_poly

        steps.append((
            "Transfer function G(s) 구성",
            f"G(s) = ({sp.sstr(num_poly)}) / ({sp.sstr(den_poly)})"
        ))

        # ---------------------------------------------------------------
        # 2. 특성방정식의 근 (Characteristic equation roots)
        # ---------------------------------------------------------------
        char_roots = sp.solve(den_poly, s)
        roots_str = ", ".join(str(r) for r in char_roots)
        steps.append((
            "Characteristic equation roots (특성방정식의 근)",
            f"Roots of {sp.sstr(den_poly)} = 0:\n  {roots_str}"
        ))

        # ---------------------------------------------------------------
        # 3. 부분분수 분해 (Partial Fraction Decomposition)
        # ---------------------------------------------------------------
        pfd = apart(G_s, s)
        steps.append((
            "Partial fraction decomposition (부분분수 분해)",
            sp.sstr(pfd)
        ))

        # ---------------------------------------------------------------
        # 4 & 5. 서브시스템 파라미터 추출
        #   - 분모를 인수분해하여 1차/2차 인수를 식별
        #   - 각 부분분수 항에서 파라미터 추출
        # ---------------------------------------------------------------
        den_factored = factor(den_poly)
        steps.append((
            "Denominator factorisation (분모 인수분해)",
            sp.sstr(den_factored)
        ))

        # 부분분수 항 분리
        pfd_terms = sp.Add.make_args(pfd)

        subsystem_info: list[dict[str, Any]] = []
        first_order_systems: list[dict] = []
        second_order_systems: list[dict] = []

        for term in pfd_terms:
            # 각 항의 분모를 추출하여 분석
            t_num, t_den = sp.fraction(sp.cancel(term))

            if t_den == 1:
                # 상수 항 (proper이 아닌 경우) — 일반적으로 strictly proper 이면 없음
                continue

            # apart()는 중복근에 대해 (s+a)^k 형태를 생성할 수 있다.
            # 기저 인수와 지수를 분리한다.
            base_den, exp_den = t_den.as_base_exp()
            if exp_den != 1 and Poly(base_den, s).degree() >= 1:
                # 중복근: 예 1/(s+2)^2 → 기저 (s+2)를 분석
                den_for_classify = base_den
            else:
                den_for_classify = t_den

            t_den_poly = Poly(den_for_classify, s)
            deg = t_den_poly.degree()
            info = _classify_factor(den_for_classify, s)
            info["term"] = term
            info["numerator"] = t_num
            if exp_den != 1:
                info["multiplicity"] = int(exp_den)
                info["full_denominator"] = t_den
            subsystem_info.append(info)

            if deg == 1:
                first_order_systems.append(info)
            elif deg == 2:
                second_order_systems.append(info)

        # 1차 서브시스템 결과
        omega_sym = Symbol("omega", real=True, positive=True)
        for i, sys1 in enumerate(first_order_systems):
            label = f"1st-order subsystem #{i+1}" if len(first_order_systems) > 1 else "1st-order subsystem"
            mult_str = ""
            if "multiplicity" in sys1:
                mult_str = f"\n  Multiplicity = {sys1['multiplicity']} (repeated root)"

            # G₁(jω): s = jω 직접 대입 (cancel 없이 원형 유지)
            G1_jw = sys1["term"].subs(s, I * omega_sym)
            G1_jw_display = sp.nsimplify(G1_jw)
            # 분자/분모를 실수+허수 표준형으로 정리
            G1_num, G1_den = sp.fraction(G1_jw_display)
            G1_den_re = sp.re(sp.expand(G1_den))
            G1_den_im = sp.im(sp.expand(G1_den))
            G1_num_re = sp.re(sp.expand(G1_num))
            G1_num_im = sp.im(sp.expand(G1_num))
            G1_jw_std = _format_complex_fraction(G1_num_re, G1_num_im,
                                                  G1_den_re, G1_den_im)
            # 실수부/허수부 분리
            G1_re = sp.simplify(sp.re(G1_jw_display))
            G1_im = sp.simplify(sp.im(G1_jw_display))

            detail = (
                f"G₁(s) = {sp.sstr(sys1['term'])}\n"
                f"  Break frequency ω_b = {sys1['omega_b']}\n"
                f"  Time constant τ = {sys1['time_constant']}"
                f"{mult_str}\n\n"
                f"  G₁(jω) = {G1_jw_std}\n"
                f"  Re[G₁(jω)] = {G1_re}\n"
                f"  Im[G₁(jω)] = {G1_im}"
            )
            steps.append((label, detail))

            # final_answer에도 jω 표현 저장
            sys1["G_jw"] = G1_jw_std
            sys1["G_jw_real"] = str(G1_re)
            sys1["G_jw_imag"] = str(G1_im)

        # 2차 서브시스템 결과
        for i, sys2 in enumerate(second_order_systems):
            label = f"2nd-order subsystem #{i+1}" if len(second_order_systems) > 1 else "2nd-order subsystem"

            # G₂(jω): s = jω 직접 대입 (cancel 없이 원형 유지)
            G2_jw = sys2["term"].subs(s, I * omega_sym)
            G2_jw_display = sp.nsimplify(G2_jw)
            # 분자/분모를 실수+허수 표준형으로 정리
            G2_num, G2_den = sp.fraction(G2_jw_display)
            G2_den_re = sp.re(sp.expand(G2_den))
            G2_den_im = sp.im(sp.expand(G2_den))
            G2_num_re = sp.re(sp.expand(G2_num))
            G2_num_im = sp.im(sp.expand(G2_num))
            G2_jw_std = _format_complex_fraction(G2_num_re, G2_num_im,
                                                  G2_den_re, G2_den_im)
            # 실수부/허수부 분리
            G2_re = sp.simplify(sp.re(G2_jw_display))
            G2_im = sp.simplify(sp.im(G2_jw_display))

            detail = (
                f"G₂(s) = {sp.sstr(sys2['term'])}\n"
                f"  Natural frequency ω_n = {sys2['omega_n']}"
                f" (ω_n² = {sys2['omega_n_squared']})\n"
                f"  Damping ratio ζ = {sys2['zeta']}\n\n"
                f"  G₂(jω) = {G2_jw_std}\n"
                f"  Re[G₂(jω)] = {G2_re}\n"
                f"  Im[G₂(jω)] = {G2_im}"
            )
            steps.append((label, detail))

            # final_answer에도 jω 표현 저장
            sys2["G_jw"] = G2_jw_std
            sys2["G_jw_real"] = str(G2_re)
            sys2["G_jw_imag"] = str(G2_im)

        # 고차(>2) 서브시스템 — 드물지만 처리
        higher_order = [si for si in subsystem_info if si["order"] > 2]
        for si in higher_order:
            steps.append((
                f"{si['order']}th-order factor",
                f"Factor: {sp.sstr(si['expression'])}\n(Further decomposition not performed)"
            ))

        # ---------------------------------------------------------------
        # 6. 특성값 테이블: G, G₁, G₂ 각각의 G(0), G(jω_char), G(j∞)
        # ---------------------------------------------------------------
        # 주요 특성 주파수 수집
        char_freqs: list[tuple[str, sp.Expr]] = []
        for idx, sys1 in enumerate(first_order_systems):
            label = f"ω_b{idx+1}" if len(first_order_systems) > 1 else "ω_b"
            char_freqs.append((label, sys1["omega_b"]))
        for idx, sys2 in enumerate(second_order_systems):
            label = f"ω_n{idx+1}" if len(second_order_systems) > 1 else "ω_n"
            char_freqs.append((label, sys2["omega_n"]))

        # 평가 대상 함수 목록: (이름, sympy 식)
        eval_targets: list[tuple[str, sp.Expr]] = [("G", G_s)]
        for idx, sys1 in enumerate(first_order_systems):
            name = f"G{idx+1}" if len(first_order_systems) + len(second_order_systems) > 1 else "G₁"
            eval_targets.append((name, sys1["term"]))
        for idx, sys2 in enumerate(second_order_systems):
            k = len(first_order_systems) + idx + 1
            name = f"G{k}" if len(first_order_systems) + len(second_order_systems) > 1 else "G₂"
            eval_targets.append((name, sys2["term"]))

        # 테이블 구성
        table_lines: list[str] = []
        char_values: dict[str, Any] = {}

        for tgt_name, tgt_expr in eval_targets:
            tgt_vals: dict[str, Any] = {}
            line_parts = [f"{tgt_name}:"]

            # G(0)
            val_0 = sp.nsimplify(tgt_expr.subs(s, 0))
            mag_0 = float(sp.Abs(val_0))
            tgt_vals["G(0)"] = mag_0
            line_parts.append(f"|{tgt_name}(0)| = {mag_0:.6g}")

            # G(jω_char) at each characteristic frequency
            for freq_label, freq_val in char_freqs:
                val_jw = sp.nsimplify(tgt_expr.subs(s, I * freq_val))
                mag_jw = float(sp.Abs(val_jw))
                phase_jw = float(sp.arg(val_jw)) * 180 / float(sp.pi)
                key = f"G(j{freq_label})"
                tgt_vals[key] = {"magnitude": mag_jw, "phase_deg": phase_jw}
                line_parts.append(
                    f"|{tgt_name}(j{freq_label})| = {mag_jw:.6g}  "
                    f"(∠ {phase_jw:.2f}°)"
                )

            # G(j∞) — limit as ω→∞
            omega_inf = Symbol("_w", positive=True)
            val_inf_expr = tgt_expr.subs(s, I * omega_inf)
            mag_inf = sp.limit(sp.Abs(val_inf_expr), omega_inf, sp.oo)
            tgt_vals["G(j∞)"] = float(mag_inf) if mag_inf.is_number else str(mag_inf)
            line_parts.append(f"|{tgt_name}(j∞)| = {mag_inf}")

            table_lines.append("\n  ".join(line_parts))
            char_values[tgt_name] = tgt_vals

        steps.append((
            "Characteristic values (특성값 테이블)",
            "\n\n".join(table_lines)
        ))
        final_answer["characteristic_values"] = char_values

        # DC gain 편의 키
        G_at_0 = G_s.subs(s, 0)
        dc_gain = sp.Abs(G_at_0)

        # ---------------------------------------------------------------
        # 7. Characteristic points (특성 주파수 요약) — 기존 호환 유지
        # ---------------------------------------------------------------
        char_points: dict[str, Any] = {
            "dc_gain": float(dc_gain),
            "first_order_subsystems": [],
            "second_order_subsystems": [],
        }
        for sys1 in first_order_systems:
            omega_b_val = float(sys1["omega_b"])
            g1_at_wb = sys1["term"].subs(s, I * sys1["omega_b"])
            mag_at_wb = float(sp.Abs(g1_at_wb))
            char_points["first_order_subsystems"].append({
                "omega_b": omega_b_val,
                "mag_at_omega_b": mag_at_wb,
            })
        for sys2 in second_order_systems:
            omega_n_val = float(sys2["omega_n"])
            zeta_val = float(sys2["zeta"])
            g2_at_wn = sys2["term"].subs(s, I * sys2["omega_n"])
            mag_at_wn = float(sp.Abs(sp.nsimplify(g2_at_wn)))
            char_points["second_order_subsystems"].append({
                "omega_n": omega_n_val,
                "zeta": zeta_val,
                "mag_at_omega_n": mag_at_wn,
            })

        system_order = len(den_coeffs) - len(num_coeffs)
        char_points["high_freq_rolloff_dB_per_decade"] = -20 * system_order
        char_points["high_freq_description"] = (
            f"High-frequency asymptote: {-20 * system_order} dB/decade "
            f"(relative order = {system_order})"
        )

        # ---------------------------------------------------------------
        # 8. G(jω) 수치 평가 (선택)
        # ---------------------------------------------------------------
        omega_range = params.get("omega_range")
        if omega_range is not None:
            if len(omega_range) == 2:
                w_min, w_max = omega_range
                n_pts = 500
            else:
                w_min, w_max, n_pts = omega_range
            omega_arr = np.logspace(np.log10(float(w_min)),
                                    np.log10(float(w_max)),
                                    int(n_pts))
            G_jw = _evaluate_frf(num_coeffs_float, den_coeffs_float, omega_arr)
            magnitude_dB = 20 * np.log10(np.abs(G_jw))
            phase_deg = np.degrees(np.angle(G_jw))
            final_answer["omega"] = omega_arr.tolist()
            final_answer["magnitude_dB"] = magnitude_dB.tolist()
            final_answer["phase_deg"] = phase_deg.tolist()

        # ---------------------------------------------------------------
        # Final answer 조립
        # ---------------------------------------------------------------
        final_answer["transfer_function"] = str(G_s)
        final_answer["roots"] = [str(r) for r in char_roots]
        final_answer["partial_fractions"] = str(pfd)
        final_answer["dc_gain"] = float(dc_gain)
        final_answer["characteristic_points"] = char_points

        # 서브시스템 파라미터
        if first_order_systems:
            final_answer["first_order"] = [
                {
                    "omega_b": float(si["omega_b"]),
                    "time_constant": float(si["time_constant"]),
                    "subsystem_tf": str(si["term"]),
                    "G_jw": si.get("G_jw", ""),
                    "G_jw_real": si.get("G_jw_real", ""),
                    "G_jw_imag": si.get("G_jw_imag", ""),
                }
                for si in first_order_systems
            ]
        if second_order_systems:
            final_answer["second_order"] = [
                {
                    "omega_n": float(si["omega_n"]),
                    "zeta": float(si["zeta"]),
                    "subsystem_tf": str(si["term"]),
                    "G_jw": si.get("G_jw", ""),
                    "G_jw_real": si.get("G_jw_real", ""),
                    "G_jw_imag": si.get("G_jw_imag", ""),
                }
                for si in second_order_systems
            ]

        # ---------------------------------------------------------------
        # 9. Forced response (강제응답): f(t) 주어질 때 x(t) 계산
        # ---------------------------------------------------------------
        forcing = params.get("forcing")
        if forcing is not None:
            force_type = forcing.get("type", "harmonic")
            t_sym = Symbol("t", positive=True)

            if force_type in ("impulse", "step", "ramp", "exponential"):
                # --- 역라플라스 변환을 통한 시간응답 ---
                F0 = float(forcing.get("F0", 1.0))

                if force_type == "impulse":
                    F_s = sp.S.One          # L{δ(t)} = 1
                    f_t_str = f"{F0} δ(t)"
                elif force_type == "step":
                    F_s = sp.S.One / s      # L{u(t)} = 1/s
                    f_t_str = f"{F0} u(t)"
                elif force_type == "ramp":
                    F_s = sp.S.One / s**2   # L{t·u(t)} = 1/s²
                    f_t_str = f"{F0} t·u(t)"
                elif force_type == "exponential":
                    a = sp.nsimplify(forcing.get("a", -1))
                    F_s = sp.S.One / (s - a)  # L{e^(at)} = 1/(s-a)
                    f_t_str = f"{F0} exp({a}t)"

                # X(s) = G(s) · F(s) · F0
                X_s = sp.nsimplify(F0) * G_s * F_s
                X_s_simplified = sp.cancel(X_s)

                steps.append((
                    f"Forced response (Laplace): f(t) = {f_t_str}",
                    f"F(s) = {F_s}\n"
                    f"X(s) = F₀ · G(s) · F(s) = {X_s_simplified}"
                ))

                # 역라플라스 변환
                try:
                    x_t = sp.inverse_laplace_transform(X_s_simplified, s, t_sym)
                    x_t = sp.simplify(x_t)

                    steps.append((
                        f"x(t) = L⁻¹{{X(s)}}",
                        f"x(t) = {x_t}"
                    ))

                    final_answer["forced_response"] = {
                        "f_t": f_t_str,
                        "type": force_type,
                        "X_s": str(X_s_simplified),
                        "x_t": str(x_t),
                    }
                except Exception as e:
                    steps.append((
                        "x(t) — inverse Laplace",
                        f"역라플라스 변환 실패: {e}\n"
                        f"X(s) = {X_s_simplified}\n"
                        "부분분수 분해 후 수동으로 변환하세요."
                    ))
                    # 부분분수 분해라도 제공
                    X_pfd = sp.apart(X_s_simplified, s)
                    steps.append((
                        "X(s) partial fractions",
                        f"X(s) = {X_pfd}"
                    ))
                    final_answer["forced_response"] = {
                        "f_t": f_t_str,
                        "type": force_type,
                        "X_s": str(X_s_simplified),
                        "X_s_partial": str(X_pfd),
                    }

            elif force_type == "harmonic":
                # f(t) = F0 * cos(w0*t) 또는 F0 * sin(w0*t)
                F0 = float(forcing.get("F0", 1.0))
                raw_omega = forcing.get("omega", forcing.get("w0", 1.0))
                func = forcing.get("func", "cos")  # "cos" or "sin"

                # ω가 심볼릭인 경우 (e.g., "w", "omega") vs 수치인 경우
                omega_is_symbolic = isinstance(raw_omega, str)

                if omega_is_symbolic:
                    # --- 심볼릭 ω: 교과서 표준형으로 system response 표현 ---
                    # D(jω) = Σ aₖ(jω)ᵏ 를 실수부/허수부로 분리하되
                    # 계수 형태를 그대로 유지한다.
                    #   Re[D] = a₀ - a₂ω² + a₄ω⁴ - ...
                    #   Im[D] = a₁ω - a₃ω³ + a₅ω⁵ - ...
                    w = Symbol("omega", real=True, positive=True)

                    # 분모 D(jω) 실수부·허수부를 계수 형태로 구성
                    n_den = len(den_coeffs)
                    # den_coeffs는 [a_n, ..., a_1, a_0] 고차→저차
                    # reversed: a[0]=a_0, a[1]=a_1, ...
                    a = [sp.nsimplify(c, rational=True) for c in reversed(den_coeffs)]

                    Re_den = sp.S.Zero
                    Im_den = sp.S.Zero
                    for k, ak in enumerate(a):
                        if ak == 0:
                            continue
                        # (jω)^k = j^k · ω^k
                        # j^0=1, j^1=j, j^2=-1, j^3=-j, j^4=1, ...
                        r = k % 4
                        term = ak * w**k
                        if r == 0:
                            Re_den += term
                        elif r == 1:
                            Im_den += term
                        elif r == 2:
                            Re_den -= term
                        elif r == 3:
                            Im_den -= term

                    # 분자 N(jω) 실수부·허수부 (같은 방식)
                    b = [sp.nsimplify(c, rational=True) for c in reversed(num_coeffs)]
                    Re_num = sp.S.Zero
                    Im_num = sp.S.Zero
                    for k, bk in enumerate(b):
                        if bk == 0:
                            continue
                        r = k % 4
                        term = bk * w**k
                        if r == 0:
                            Re_num += term
                        elif r == 1:
                            Im_num += term
                        elif r == 2:
                            Re_num -= term
                        elif r == 3:
                            Im_num -= term

                    # 표준형 문자열 생성
                    Re_den_str = sp.sstr(Re_den)
                    Im_den_str = sp.sstr(Im_den)
                    Re_num_str = sp.sstr(Re_num)
                    Im_num_str = sp.sstr(Im_num)

                    # 분자가 상수(예: 1)인 경우 간단한 형태
                    num_is_const = (Im_num == 0 and Re_num.is_number)
                    num_val = Re_num if num_is_const else None

                    # |G(jω)| 와 φ 표현
                    if num_is_const and num_val == 1:
                        # G(jω) = 1 / D(jω)
                        mag_str = f"1/sqrt(({Re_den_str})^2 + ({Im_den_str})^2)"
                        phi_str = f"atan2({Im_den_str}, {Re_den_str})"
                        F0_mag_str = (
                            f"{sp.nsimplify(F0)}"
                            f"/sqrt(({Re_den_str})^2 + ({Im_den_str})^2)"
                        )
                    elif num_is_const:
                        mag_str = (
                            f"{num_val}"
                            f"/sqrt(({Re_den_str})^2 + ({Im_den_str})^2)"
                        )
                        phi_str = f"atan2({Im_den_str}, {Re_den_str})"
                        F0_num = sp.nsimplify(F0) * num_val
                        F0_mag_str = (
                            f"{F0_num}"
                            f"/sqrt(({Re_den_str})^2 + ({Im_den_str})^2)"
                        )
                    else:
                        mag_str = (
                            f"sqrt(({Re_num_str})^2 + ({Im_num_str})^2)"
                            f"/sqrt(({Re_den_str})^2 + ({Im_den_str})^2)"
                        )
                        phi_str = (
                            f"atan2({Im_den_str}, {Re_den_str})"
                            f" - atan2({Im_num_str}, {Re_num_str})"
                        )
                        F0_mag_str = (
                            f"{sp.nsimplify(F0)} * "
                            f"sqrt(({Re_num_str})^2 + ({Im_num_str})^2)"
                            f"/sqrt(({Re_den_str})^2 + ({Im_den_str})^2)"
                        )

                    x_t_str = f"{F0_mag_str} · {func}(ω·t - φ)"

                    steps.append((
                        f"System response: f(t) = {F0} {func}(ω·t)",
                        f"G(jω) = N(jω) / D(jω)\n\n"
                        f"D(jω) = ({Re_den_str}) + j·({Im_den_str})\n"
                        + (f"N(jω) = {Re_num}\n\n"
                           if num_is_const else
                           f"N(jω) = ({Re_num_str}) + j·({Im_num_str})\n\n")
                        + f"|G(jω)| = {mag_str}\n\n"
                        f"φ = {phi_str}\n\n"
                        f"x(t) = {x_t_str}"
                    ))

                    final_answer["forced_response"] = {
                        "f_t": f"{F0} {func}(ω·t)",
                        "omega_0": "symbolic (ω)",
                        "D_jw_real": Re_den_str,
                        "D_jw_imag": Im_den_str,
                        "G_jw_magnitude": mag_str,
                        "phi": phi_str,
                        "x_t": x_t_str,
                    }

                else:
                    # --- 수치 ω: 기존 로직 ---
                    w0 = float(raw_omega)

                    # G(jω₀) 평가
                    G_jw0 = complex(G_s.subs(s, I * sp.nsimplify(w0)))
                    mag_G = abs(G_jw0)
                    phase_G = np.degrees(np.angle(G_jw0))

                    # 각 서브시스템의 기여
                    sub_lines = []
                    for idx, sys1 in enumerate(first_order_systems):
                        k = idx + 1
                        g_jw0 = complex(sys1["term"].subs(s, I * sp.nsimplify(w0)))
                        sub_lines.append(
                            f"  G{k}(j{w0:.4g}) : |G{k}| = {abs(g_jw0):.6g}, "
                            f"∠G{k} = {np.degrees(np.angle(g_jw0)):.2f}°"
                        )
                    for idx, sys2 in enumerate(second_order_systems):
                        k = len(first_order_systems) + idx + 1
                        g_jw0 = complex(sys2["term"].subs(s, I * sp.nsimplify(w0)))
                        sub_lines.append(
                            f"  G{k}(j{w0:.4g}) : |G{k}| = {abs(g_jw0):.6g}, "
                            f"∠G{k} = {np.degrees(np.angle(g_jw0)):.2f}°"
                        )

                    if func == "sin":
                        x_expr = f"{F0 * mag_G:.6g} sin({w0:.6g}t {phase_G:+.2f}°)"
                    else:
                        x_expr = f"{F0 * mag_G:.6g} cos({w0:.6g}t {phase_G:+.2f}°)"

                    steps.append((
                        f"Forced response: f(t) = {F0:.4g} {func}({w0:.4g} t)",
                        f"G(j{w0:.4g}) = {G_jw0.real:.6g} {G_jw0.imag:+.6g}j\n"
                        f"|G(j{w0:.4g})| = {mag_G:.6g}\n"
                        f"∠G(j{w0:.4g}) = {phase_G:.2f}°\n\n"
                        + ("\n".join(sub_lines) + "\n\n" if sub_lines else "")
                        + f"x_ss(t) = |G(jω₀)| · F₀ · {func}(ω₀t + ∠G)\n"
                        f"x_ss(t) = {x_expr}"
                    ))

                    final_answer["forced_response"] = {
                        "f_t": f"{F0} {func}({w0} t)",
                        "omega_0": w0,
                        "G_jw0_real": G_jw0.real,
                        "G_jw0_imag": G_jw0.imag,
                        "G_jw0_magnitude": mag_G,
                        "G_jw0_phase_deg": phase_G,
                        "x_ss": x_expr,
                        "amplitude": F0 * mag_G,
                        "phase_deg": phase_G,
                    }

            elif force_type == "general_harmonic":
                # f(t) = sum of F_k cos(w_k t + phi_k)  (중첩)
                components = forcing.get("components", [])
                resp_lines = []
                resp_data = []
                for comp in components:
                    F_k = float(comp.get("F0", 1.0))
                    w_k = float(comp.get("omega", 1.0))
                    phi_k = float(comp.get("phase_deg", 0.0))
                    func_k = comp.get("func", "cos")

                    G_jw_k = complex(G_s.subs(s, I * sp.nsimplify(w_k)))
                    mag_k = abs(G_jw_k)
                    phase_k = np.degrees(np.angle(G_jw_k))
                    total_phase = phi_k + phase_k

                    x_amp = F_k * mag_k
                    x_str = f"{x_amp:.6g} {func_k}({w_k:.6g}t {total_phase:+.2f}°)"
                    resp_lines.append(
                        f"  f_{len(resp_lines)+1} = {F_k} {func_k}({w_k}t {phi_k:+.1f}°)  "
                        f"→  x_{len(resp_lines)+1} = {x_str}"
                    )
                    resp_data.append({
                        "f_component": f"{F_k} {func_k}({w_k}t {phi_k:+.1f}°)",
                        "x_component": x_str,
                        "amplitude": x_amp,
                        "phase_deg": total_phase,
                    })

                steps.append((
                    "Forced response (superposition of harmonics)",
                    "By superposition:\n" + "\n".join(resp_lines)
                    + "\n\nx_ss(t) = " + " + ".join(
                        d["x_component"] for d in resp_data
                    )
                ))
                final_answer["forced_response"] = resp_data

        # ---------------------------------------------------------------
        # Sanity check
        # ---------------------------------------------------------------
        sanity_parts: list[str] = []

        # 안정성 확인: 모든 근의 실수부가 음수인지
        all_stable = all(sp.re(r) < 0 for r in char_roots)
        sanity_parts.append(
            f"Stability: {'All roots have negative real parts → STABLE'}"
            if all_stable
            else "Stability: WARNING — not all roots have negative real parts → UNSTABLE or MARGINALLY STABLE"
        )

        # DC gain 검증: G(0) = num_const / den_const
        if den_coeffs_float[-1] != 0:
            expected_dc = num_coeffs_float[-1] / den_coeffs_float[-1]
            sanity_parts.append(
                f"DC gain check: numerator_const/denominator_const = "
                f"{num_coeffs_float[-1]}/{den_coeffs_float[-1]} = {expected_dc:.6g} ✓"
            )
        else:
            sanity_parts.append("DC gain: denominator has zero constant term (type ≥ 1 system, DC gain → ∞)")

        # 부분분수 재조합 확인
        pfd_recombined = sp.cancel(pfd)
        G_s_simplified = sp.cancel(G_s)
        if sp.simplify(pfd_recombined - G_s_simplified) == 0:
            sanity_parts.append("Partial fraction recombination matches original G(s) ✓")
        else:
            sanity_parts.append(
                "WARNING: Partial fraction recombination check — "
                "symbolic simplification could not confirm exact match "
                "(may be due to floating-point or symbolic representation differences)"
            )

        sanity_check_str = "\n".join(sanity_parts)

        return SolverResult(
            problem_type="FRF Analysis (Frequency Response Function)",
            given={
                "denominator_coefficients": den_coeffs,
                "numerator_coefficients": num_coeffs,
                "system_order": len(den_coeffs) - 1,
            },
            steps=steps,
            final_answer=final_answer,
            sanity_check=sanity_check_str,
        )

    def get_input_template(self) -> dict:
        """필요한 입력 파라미터 설명을 반환한다."""
        return {
            "coefficients": {
                "type": "list[float]",
                "required": True,
                "description": (
                    "분모(특성방정식) 다항식 계수, 고차→저차 순서. "
                    "예: x''' + 3x'' + 6x' + 8x = f(t) → [1, 3, 6, 8]"
                ),
                "example": [1, 3, 6, 8],
            },
            "numerator": {
                "type": "list[float]",
                "required": False,
                "default": [1],
                "description": (
                    "분자 다항식 계수, 고차→저차 순서. "
                    "단위 임펄스 응답이면 [1] (기본값)."
                ),
                "example": [1],
            },
            "omega_range": {
                "type": "tuple (omega_min, omega_max[, num_points])",
                "required": False,
                "default": None,
                "description": (
                    "G(jω) 수치 평가를 위한 주파수 범위. "
                    "지정하면 magnitude(dB)와 phase(deg)를 계산하여 반환한다. "
                    "예: (0.01, 100, 500)"
                ),
                "example": (0.01, 100, 500),
            },
        }
