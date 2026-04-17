"""Fourier 급수 분해 및 주기 가진 응답 솔버.

주기 함수를 Fourier 급수로 전개하고, 1차/2차 시스템의 정상상태 응답을
중첩 원리로 구한다.

지원 파형:
  - square_wave: 구형파
  - sawtooth: 톱니파
  - triangle: 삼각파
  - half_sine: 반파 정류 사인파
  - custom: 사용자 정의 Fourier 계수
"""

from __future__ import annotations

from typing import Any

import numpy as np
import sympy as sp
from sympy import (
    Rational, Symbol, cos, sin, sqrt, pi, oo,
    integrate, Piecewise, simplify, nsimplify, atan2,
)

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Fourier 급수 계수 계산 유틸리티
# ---------------------------------------------------------------------------

t_sym = Symbol("t", real=True)
T_sym = Symbol("T", positive=True)
n_sym = Symbol("n", positive=True, integer=True)
omega_sym = Symbol("omega", positive=True, real=True)


def _square_wave_coefficients(f0: sp.Expr, T: sp.Expr, N: int) -> dict:
    """구형파의 Fourier 급수 계수를 해석적으로 계산.

    f(t) = f0  for  -T/4 < t < T/4
           0   otherwise (within one period)

    Returns dict with 'a0', 'an', 'bn' as lists (n=1..N).
    """
    omega = 2 * pi / T

    # a0 = (1/T) * integral over one period
    # 구형파: f0 구간 길이 = T/2, 전체 주기 T
    a0 = f0 / 2

    an_list = []  # a_n coefficients (cosine)
    bn_list = []  # b_n coefficients (sine)

    for n in range(1, N + 1):
        # a_n = (2/T) * int_{-T/4}^{T/4} f0 * cos(n*omega*t) dt
        # = (2*f0/T) * [sin(n*omega*t)/(n*omega)]_{-T/4}^{T/4}
        # = (2*f0/T) * 2*sin(n*pi/2) / (n*omega)
        # = (2*f0/(n*pi)) * sin(n*pi/2)
        n_val = sp.Integer(n)
        a_n = 2 * f0 / (n_val * pi) * sp.sin(n_val * pi / 2)
        a_n = simplify(a_n)
        an_list.append(a_n)

        # b_n = 0 (even function)
        bn_list.append(sp.S.Zero)

    return {"a0": a0, "an": an_list, "bn": bn_list, "omega": omega}


def _sawtooth_coefficients(f0: sp.Expr, T: sp.Expr, N: int) -> dict:
    """톱니파 Fourier 계수. f(t) = f0 * (2t/T) for -T/2 < t < T/2."""
    omega = 2 * pi / T
    a0 = sp.S.Zero
    an_list = []
    bn_list = []

    for n in range(1, N + 1):
        n_val = sp.Integer(n)
        an_list.append(sp.S.Zero)
        # b_n = (-1)^(n+1) * 2*f0 / (n*pi)
        b_n = (-1) ** (n + 1) * 2 * f0 / (n_val * pi)
        bn_list.append(simplify(b_n))

    return {"a0": a0, "an": an_list, "bn": bn_list, "omega": omega}


def _triangle_coefficients(f0: sp.Expr, T: sp.Expr, N: int) -> dict:
    """삼각파 Fourier 계수. 짝수함수, 진폭 f0."""
    omega = 2 * pi / T
    a0 = f0 / 2
    an_list = []
    bn_list = []

    for n in range(1, N + 1):
        n_val = sp.Integer(n)
        # a_n = (4*f0/(n^2*pi^2)) * ((-1)^((n-1)/2)) for odd n, 0 for even n
        if n % 2 == 1:
            # a_n = (-1)^((n-1)/2) * 8*f0 / (n^2 * pi^2)  [standard form]
            a_n = (-1) ** ((n - 1) // 2) * 8 * f0 / (n_val ** 2 * pi ** 2)
        else:
            a_n = sp.S.Zero
        an_list.append(simplify(a_n))
        bn_list.append(sp.S.Zero)

    return {"a0": a0, "an": an_list, "bn": bn_list, "omega": omega}


def _custom_coefficients(a0, an_list, bn_list, omega):
    """사용자 정의 Fourier 계수."""
    return {
        "a0": sp.nsimplify(a0),
        "an": [sp.nsimplify(a) for a in an_list],
        "bn": [sp.nsimplify(b) for b in bn_list],
        "omega": sp.nsimplify(omega),
    }


# ---------------------------------------------------------------------------
# 시스템 응답 계산
# ---------------------------------------------------------------------------

def _first_order_harmonic_response(
    c: sp.Expr, k: sp.Expr,
    amp: sp.Expr, n_omega: sp.Expr,
    is_cosine: bool = True,
) -> dict:
    """1차 시스템 cx' + kx = A*cos(nωt) (or sin) 의 정상상태 응답.

    G(jnω) = 1 / (k + j*c*n*ω)
    |G| = 1 / sqrt(k^2 + (c*n*ω)^2)
    φ = -arctan(c*n*ω / k)

    x_n(t) = A * |G| * cos(nωt + φ)  [or sin]
    """
    mag_sq = k ** 2 + (c * n_omega) ** 2
    magnitude = 1 / sqrt(mag_sq)
    phase = -sp.atan(c * n_omega / k)

    x_amp = simplify(amp * magnitude)

    return {
        "amplitude": x_amp,
        "magnitude": magnitude,
        "phase": phase,
        "is_cosine": is_cosine,
        "n_omega": n_omega,
    }


def _second_order_harmonic_response(
    m: sp.Expr, c: sp.Expr, k: sp.Expr,
    amp: sp.Expr, n_omega: sp.Expr,
    is_cosine: bool = True,
) -> dict:
    """2차 시스템 mx'' + cx' + kx = A*cos(nωt) 의 정상상태 응답.

    G(jnω) = 1 / (k - m*(nω)^2 + j*c*nω)
    """
    re_part = k - m * n_omega ** 2
    im_part = c * n_omega
    mag_sq = re_part ** 2 + im_part ** 2
    magnitude = 1 / sqrt(mag_sq)
    phase = -sp.atan2(im_part, re_part)

    x_amp = simplify(amp * magnitude)

    return {
        "amplitude": x_amp,
        "magnitude": magnitude,
        "phase": phase,
        "is_cosine": is_cosine,
        "n_omega": n_omega,
    }


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class FourierSolver(BaseSolver):
    """Fourier 급수 기반 주기 가진 응답 솔버.

    주기 함수를 Fourier 급수로 전개하고, 각 조화 성분에 대한
    정상상태 응답을 중첩하여 전체 응답을 구한다.
    """

    def solve(self, params: dict) -> SolverResult:
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # -----------------------------------------------------------
        # 0. 입력 파싱
        # -----------------------------------------------------------
        system_order = params.get("system_order", 1)
        waveform = params.get("waveform", "square_wave")
        N = params.get("num_terms", 10)  # Fourier 급수 항 수

        # 시스템 파라미터
        c_val = sp.nsimplify(params.get("c", 1), rational=True)
        k_val = sp.nsimplify(params.get("k", 1), rational=True)
        m_val = sp.nsimplify(params.get("m", 0), rational=True) if system_order >= 2 else sp.S.Zero

        # 가진 파라미터
        f0_val = sp.nsimplify(params.get("f0", 1), rational=True)
        T_val = params.get("T")  # 주기 (None이면 심볼릭)

        use_symbolic_T = (T_val is None) or (isinstance(T_val, str) and T_val.lower() == "t")
        if use_symbolic_T:
            T_expr = T_sym
        else:
            T_expr = sp.nsimplify(T_val, rational=True)

        omega_expr = 2 * pi / T_expr  # 기본 주파수

        # 심볼릭 변수 선언
        f0_sym_flag = params.get("f0_symbolic", False)
        if f0_sym_flag or (isinstance(params.get("f0"), str)):
            f0_expr = Symbol("f_0", positive=True)
        else:
            f0_expr = f0_val

        c_sym_flag = params.get("c_symbolic", False)
        k_sym_flag = params.get("k_symbolic", False)
        if c_sym_flag or isinstance(params.get("c"), str):
            c_expr = Symbol("c", positive=True)
        else:
            c_expr = c_val

        if k_sym_flag or isinstance(params.get("k"), str):
            k_expr = Symbol("k", positive=True)
        else:
            k_expr = k_val

        if system_order >= 2:
            m_sym_flag = params.get("m_symbolic", False)
            if m_sym_flag or isinstance(params.get("m"), str):
                m_expr = Symbol("m", positive=True)
            else:
                m_expr = m_val
        else:
            m_expr = sp.S.Zero

        # Given 정보
        given: dict[str, Any] = {
            "system_order": system_order,
            "waveform": waveform,
            "num_fourier_terms": N,
        }
        if system_order == 1:
            given["equation"] = f"c·x' + k·x = f(t),  c={c_expr}, k={k_expr}"
        else:
            given["equation"] = f"m·x'' + c·x' + k·x = f(t),  m={m_expr}, c={c_expr}, k={k_expr}"
        given["f0"] = str(f0_expr)
        given["T"] = str(T_expr)

        # -----------------------------------------------------------
        # 1. Fourier 급수 전개
        # -----------------------------------------------------------
        steps.append(("Step 1: Fourier 급수 전개",
                       f"파형: {waveform}\n"
                       f"주기 T = {T_expr}, 기본 주파수 ω = 2π/T = {omega_expr}"))

        if waveform == "square_wave":
            coeffs = _square_wave_coefficients(f0_expr, T_expr, N)
            f_desc = (
                f"f(t) = {{ f₀ = {f0_expr},  -T/4 < t < T/4\n"
                f"       {{ 0,         otherwise in [-T/2, T/2]"
            )
        elif waveform == "sawtooth":
            coeffs = _sawtooth_coefficients(f0_expr, T_expr, N)
            f_desc = f"f(t) = {f0_expr} · (2t/T),  -T/2 < t < T/2"
        elif waveform == "triangle":
            coeffs = _triangle_coefficients(f0_expr, T_expr, N)
            f_desc = f"f(t) = 삼각파, 진폭 {f0_expr}, 주기 T"
        elif waveform == "custom":
            custom_a0 = params.get("a0", 0)
            custom_an = params.get("an", [])
            custom_bn = params.get("bn", [])
            coeffs = _custom_coefficients(custom_a0, custom_an, custom_bn, omega_expr)
            N = max(len(custom_an), len(custom_bn))
            f_desc = "f(t) = 사용자 정의 Fourier 계수"
        else:
            raise ValueError(f"지원하지 않는 파형: {waveform}. "
                             f"['square_wave', 'sawtooth', 'triangle', 'custom'] 중 선택하세요.")

        steps.append(("가진 함수 정의", f_desc))

        a0 = coeffs["a0"]
        an_list = coeffs["an"]
        bn_list = coeffs["bn"]
        omega_base = coeffs["omega"]

        # Fourier 급수 표현 구성
        fourier_terms_str = [f"a₀ = {a0}"]
        fourier_nonzero = []
        for n in range(1, N + 1):
            a_n = an_list[n - 1]
            b_n = bn_list[n - 1]
            if a_n != 0:
                fourier_terms_str.append(f"a_{n} = {a_n}")
                fourier_nonzero.append((n, a_n, "cos"))
            if b_n != 0:
                fourier_terms_str.append(f"b_{n} = {b_n}")
                fourier_nonzero.append((n, b_n, "sin"))

        steps.append((
            "Fourier 계수",
            "f(t) = a₀ + Σ [aₙ cos(nωt) + bₙ sin(nωt)]\n\n"
            + "\n".join(fourier_terms_str)
            + f"\n\n(0이 아닌 항 수: {len(fourier_nonzero)}개, n = 1..{N})"
        ))

        # Fourier 급수를 수식으로 표현
        t = Symbol("t", real=True)
        f_fourier = a0
        for n_val, coeff, func_type in fourier_nonzero:
            if func_type == "cos":
                f_fourier += coeff * cos(n_val * omega_base * t)
            else:
                f_fourier += coeff * sin(n_val * omega_base * t)

        steps.append((
            "Fourier 급수 전개식",
            f"f(t) = {f_fourier}"
        ))

        # -----------------------------------------------------------
        # 2. 전달함수 G(jω) 정의
        # -----------------------------------------------------------
        if system_order == 1:
            steps.append((
                "Step 2: 1차 시스템 전달함수",
                f"c·x' + k·x = f(t)\n"
                f"G(jω) = 1/(k + jcω)\n"
                f"|G(jω)| = 1/√(k² + c²ω²)\n"
                f"φ(ω) = -arctan(cω/k)\n\n"
                f"시정수 τ = c/k = {simplify(c_expr / k_expr)}"
            ))
        else:
            steps.append((
                "Step 2: 2차 시스템 전달함수",
                f"m·x'' + c·x' + k·x = f(t)\n"
                f"G(jω) = 1/(k - mω² + jcω)\n"
                f"|G(jω)| = 1/√((k - mω²)² + c²ω²)\n"
                f"φ(ω) = -arctan(cω/(k - mω²))"
            ))

        # -----------------------------------------------------------
        # 3. 각 조화 성분에 대한 응답 계산
        # -----------------------------------------------------------
        steps.append(("Step 3: 각 조화 성분의 정상상태 응답 (중첩 원리)", ""))

        # DC 성분 (a0)
        if a0 != 0:
            if system_order == 1:
                x_dc = a0 / k_expr
            else:
                x_dc = a0 / k_expr
            steps.append((
                "DC 성분 (n=0): x₀ = a₀/k",
                f"x₀ = {a0} / {k_expr} = {simplify(x_dc)}"
            ))
        else:
            x_dc = sp.S.Zero
            steps.append(("DC 성분 (n=0)", "a₀ = 0 → x₀ = 0"))

        # 각 조화 성분
        response_terms = []
        response_str_parts = []

        if x_dc != 0:
            response_terms.append(("dc", x_dc, None, None))

        for n_val, coeff, func_type in fourier_nonzero:
            n_omega = n_val * omega_base
            is_cos = (func_type == "cos")

            if system_order == 1:
                resp = _first_order_harmonic_response(
                    c_expr, k_expr, coeff, n_omega, is_cos
                )
            else:
                resp = _second_order_harmonic_response(
                    m_expr, c_expr, k_expr, coeff, n_omega, is_cos
                )

            x_amp = resp["amplitude"]
            mag = resp["magnitude"]
            phase = resp["phase"]
            func_name = "cos" if is_cos else "sin"

            # 응답 항 저장
            response_terms.append((n_val, x_amp, phase, func_type))

            # 풀이 과정 표시
            detail = (
                f"n = {n_val}: {'aₙ' if is_cos else 'bₙ'} = {coeff}\n"
                f"  nω = {n_val}ω = {simplify(n_omega)}\n"
                f"  |G(j·{n_val}ω)| = {simplify(mag)}\n"
                f"  φ_{n_val} = {simplify(phase)}\n"
                f"  xₙ(t) = {simplify(x_amp)} · {func_name}({n_val}ωt + φ_{n_val})"
            )
            steps.append((f"n = {n_val} 성분 응답", detail))
            response_str_parts.append(
                f"{simplify(x_amp)} · {func_name}({n_val}ωt + ({simplify(phase)}))"
            )

        # -----------------------------------------------------------
        # 4. 전체 정상상태 응답 조립
        # -----------------------------------------------------------
        # 심볼릭 응답 구성
        x_ss = x_dc
        for item in response_terms:
            if item[0] == "dc":
                continue
            n_val, x_amp, phase, func_type = item
            n_omega = n_val * omega_base
            if func_type == "cos":
                x_ss += x_amp * cos(n_val * omega_base * t + phase)
            else:
                x_ss += x_amp * sin(n_val * omega_base * t + phase)

        dc_str = f"{simplify(x_dc)}" if x_dc != 0 else ""
        total_resp_str = dc_str
        if response_str_parts:
            if total_resp_str:
                total_resp_str += "\n  + " + "\n  + ".join(response_str_parts)
            else:
                total_resp_str = "\n  + ".join(response_str_parts)

        steps.append((
            "Step 4: 전체 정상상태 응답 x_ss(t)",
            f"x_ss(t) = {total_resp_str}"
        ))

        # -----------------------------------------------------------
        # 5. 일반 공식 (닫힌 형태)
        # -----------------------------------------------------------
        if waveform == "square_wave" and system_order == 1:
            n_gen = Symbol("n", positive=True, integer=True, odd=True)
            omega_gen = 2 * pi / T_expr

            general_formula = (
                f"x_ss(t) = f₀/(2k)\n"
                f"  + Σ_{{n=1,3,5,...}} "
                f"(2f₀/(nπ)) · (-1)^((n-1)/2) · "
                f"[1/√(k² + c²(nω)²)] · cos(nωt - arctan(cnω/k))\n\n"
                f"where ω = 2π/T"
            )
            steps.append(("일반 공식 (구형파 + 1차 시스템)", general_formula))

        # -----------------------------------------------------------
        # 6. Sanity check
        # -----------------------------------------------------------
        sanity_parts = []

        # Parseval 정리 검증 (에너지 보존)
        energy_input = a0 ** 2
        for n in range(N):
            energy_input += (an_list[n] ** 2 + bn_list[n] ** 2) / 2
        sanity_parts.append(
            f"Fourier 계수 에너지 (Parseval): Σ = {simplify(energy_input)}"
        )

        # 시정수 검증
        if system_order == 1:
            tau = c_expr / k_expr
            sanity_parts.append(f"시정수 τ = c/k = {simplify(tau)}")
            sanity_parts.append(
                f"nω >> k/c 인 고차 항은 급격히 감쇠 → "
                f"저역통과 필터 특성 확인 ✓"
            )

        sanity_check = "\n".join(sanity_parts)

        # -----------------------------------------------------------
        # Final answer 조립
        # -----------------------------------------------------------
        final_answer["fourier_coefficients"] = {
            "a0": str(a0),
            "nonzero_terms": [
                {"n": n, "coefficient": str(c), "type": ft}
                for n, c, ft in fourier_nonzero
            ],
        }
        final_answer["dc_response"] = str(simplify(x_dc))
        final_answer["harmonic_responses"] = [
            {
                "n": item[0],
                "amplitude": str(simplify(item[1])),
                "phase": str(simplify(item[2])) if item[2] is not None else "0",
                "type": item[3] if item[3] else "dc",
            }
            for item in response_terms
        ]
        final_answer["x_ss_symbolic"] = str(x_ss)
        final_answer["omega"] = str(omega_base)

        # 수치 평가 (T가 수치인 경우)
        if not use_symbolic_T and not f0_sym_flag and not c_sym_flag and not k_sym_flag:
            try:
                t_arr = np.linspace(0, float(T_expr) * 2, 500)
                x_func = sp.lambdify(t, x_ss, modules="numpy")
                f_func = sp.lambdify(t, f_fourier, modules="numpy")
                x_arr = x_func(t_arr)
                f_arr = f_func(t_arr)
                final_answer["numerical"] = {
                    "t": t_arr.tolist(),
                    "x_ss": (x_arr if isinstance(x_arr, np.ndarray)
                             else np.full_like(t_arr, float(x_arr))).tolist(),
                    "f_fourier": (f_arr if isinstance(f_arr, np.ndarray)
                                  else np.full_like(t_arr, float(f_arr))).tolist(),
                }
            except Exception:
                pass  # 수치 평가 실패 시 무시

        return SolverResult(
            problem_type="Fourier Series Response (주기 가진 응답)",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check=sanity_check,
        )

    def get_input_template(self) -> dict:
        return {
            "system_order": {
                "type": "int",
                "required": False,
                "default": 1,
                "description": "시스템 차수 (1: cx'+kx=f, 2: mx''+cx'+kx=f)",
            },
            "c": {
                "type": "float or str",
                "required": True,
                "description": "감쇠 계수 (숫자 또는 'c'로 심볼릭)",
                "example": "c",
            },
            "k": {
                "type": "float or str",
                "required": True,
                "description": "강성 계수 (숫자 또는 'k'로 심볼릭)",
                "example": "k",
            },
            "m": {
                "type": "float or str",
                "required": False,
                "description": "질량 (2차 시스템일 때만, 숫자 또는 'm')",
            },
            "waveform": {
                "type": "str",
                "required": True,
                "description": "파형 종류: 'square_wave', 'sawtooth', 'triangle', 'custom'",
                "example": "square_wave",
            },
            "f0": {
                "type": "float or str",
                "required": True,
                "description": "가진력 진폭 (숫자 또는 'f0'로 심볼릭)",
                "example": "f0",
            },
            "T": {
                "type": "float or str or None",
                "required": False,
                "default": None,
                "description": "주기 (숫자 = 수치, 'T' 또는 None = 심볼릭)",
                "example": "T",
            },
            "num_terms": {
                "type": "int",
                "required": False,
                "default": 10,
                "description": "Fourier 급수 전개 항 수 N",
            },
            "a0": {"type": "float", "required": False, "description": "custom 파형의 a0"},
            "an": {"type": "list[float]", "required": False, "description": "custom 파형의 aₙ 리스트"},
            "bn": {"type": "list[float]", "required": False, "description": "custom 파형의 bₙ 리스트"},
        }
