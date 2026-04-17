"""ME551 Midterm Solver v1.0 -- fully interactive CLI.

Wires together all solver engines and templates into a polished
exam-day command-line interface.
"""

from __future__ import annotations

import traceback
from datetime import datetime
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Solver imports (lazy-ish: all at top for clarity)
# ---------------------------------------------------------------------------
from ..core.base import SolverResult
from ..core.report import ReportEngine

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VERSION = "2.0"
_REPORTS_DIR = Path(__file__).resolve().parent.parent / "reports"

BANNER = rf"""
=============================================
   ME551 Midterm Solver v{_VERSION}
   Linear Vibration Engineering
=============================================
"""

MENU = """
------ Main Menu ----------------------------
 a. 자동 분류 (문제 입력 -> 유형 판별 -> 풀이)
 ---
 1. T/F 개념 판별
 2. FRF 분해 (전달함수)
 3. State-Space / 전이행렬 (상태방정식)
 4. Inverse Laplace Transform (역 라플라스)
 5. Lagrange -> 평형점 -> 선형화 -> 안정성
 6. N-DOF Modal Analysis
 7. Proportional Damping Response
 8. 확장 (Gyroscopic / Nonconservative)
 9. Base Excitation (지반 가진)
10. Fourier 급수 (주기 가진 응답)
11. Convolution Integral (Duhamel 적분)
 ---
12. 빈출 템플릿: Rotating Hoop
13. 빈출 템플릿: Rotating Triangle
14. 빈출 템플릿: 2-DOF Chain System
 ---
 h. 도움말 (Help)
 0. 종료
---------------------------------------------
"""

HELP_TEXT = """
=============================================
            도움말 / Help
=============================================

[1] T/F 개념 판별
    강의 내용 기반 True/False 판별.
    키워드나 문장을 입력하면 관련 규칙을
    검색하여 판정 결과와 이유를 보여줍니다.

[2] FRF 분해 (전달함수)
    ODE를 입력하면 전달함수 G(s)를 구성하고,
    부분분수 분해를 수행하여 각 서브시스템의
    ω_b, ω_n, ζ 를 추출합니다.

[3] State-Space / 전이행렬 (상태방정식)
    두 가지 모드 지원:
    a) A, B, x0 직접 입력 → ẋ = Ax + Bu
    b) M, C, K → 2차 시스템을 상태공간으로 변환
    전이행렬 Φ(t)=e^(At) 계산 (고유값분해 + Laplace),
    고유값, 안정성, 자유/강제 응답을 계산합니다.

[4] Inverse Laplace Transform (역 라플라스)
    F(s) = N(s)/D(s) 형태의 함수를 입력하면
    극점 분석, 부분분수 분해, 역 라플라스 변환을
    단계별로 수행하여 f(t)를 도출합니다.
    입력: 계수 리스트 또는 수식 문자열

[5] Lagrange -> 평형점 -> 선형화 -> 안정성
    운동에너지 T, 위치에너지 V 식을 입력하면
    Lagrange 방정식으로 EOM을 유도하고,
    평형점, 선형화, 안정성 분석까지 수행합니다.

[6] N-DOF Modal Analysis
    M, K 행렬을 입력하면 고유값/고유벡터,
    모드형상, 초기조건 응답을 계산합니다.

[7] Proportional Damping Response
    M, K, C (또는 alpha, beta) 행렬을 입력하면
    비례감쇠 모달 분석과 감쇠 응답을 계산합니다.

[a] 자동 분류
    문제 텍스트를 입력하면 유형을 자동 판별하여
    해당 솔버를 실행합니다.

[8] 확장 (Gyroscopic / Nonconservative)
    M, K, G, H, C 행렬을 입력하면
    고유값, 안정성, 모드를 분석합니다.

[9] Base Excitation (지반 가진)
    mx'' + cx' + kx = cy' + ky, y = A sin(ωt)
    전달률, 상대변위, 정상상태 응답, 시간응답을
    계산합니다. SDOF / MDOF 모두 지원.

[10] Rotating Hoop  (빈출 템플릿)
[11] Rotating Triangle  (빈출 템플릿)
[12] 2-DOF Chain System  (빈출 템플릿)

=============================================
"""

# ---------------------------------------------------------------------------
# Input parsing helpers
# ---------------------------------------------------------------------------


def parse_matrix(input_str: str) -> list[list[float]]:
    """Parse a matrix string like '9,0;0,9' into [[9,0],[0,9]].

    Rows are separated by ';', columns by ','.
    """
    input_str = input_str.strip()
    if not input_str:
        raise ValueError("빈 입력입니다.")
    rows = input_str.split(";")
    matrix = []
    for row in rows:
        row = row.strip()
        if not row:
            continue
        matrix.append([float(x.strip()) for x in row.split(",")])
    # Validate rectangular
    ncols = len(matrix[0])
    for i, row in enumerate(matrix):
        if len(row) != ncols:
            raise ValueError(f"행 {i+1}의 열 수({len(row)})가 첫 행({ncols})과 다릅니다.")
    return matrix


def parse_vector(input_str: str) -> list[float]:
    """Parse a vector string like '1,2,3' into [1.0, 2.0, 3.0]."""
    input_str = input_str.strip()
    if not input_str:
        raise ValueError("빈 입력입니다.")
    return [float(x.strip()) for x in input_str.split(",")]


def parse_float_or_symbol(input_str: str) -> float | None:
    """Return float if numeric, None if 'symbolic' or empty."""
    s = input_str.strip().lower()
    if s in ("", "s", "sym", "symbolic", "none"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parse_forcing_string(force_str: str) -> dict[str, Any]:
    """Parse a forcing function string into a forcing spec dict.

    Supported formats:
      Special signals:
        "impulse", "delta"          → impulse (Dirac delta)
        "step", "u(t)", "u"        → unit step
        "3step", "5u(t)"           → scaled step
        "ramp", "t"                → ramp t*u(t)
        "3ramp", "2t"              → scaled ramp
        "exp(-2t)", "e^(-3t)"      → exponential decay
        "5exp(-2t)"                → scaled exponential
      Harmonic:
        "cos", "3cos2", "5sin(3)"  → single harmonic
        "2cos(1)+3sin(4)"          → superposition
    """
    import re

    force_str = force_str.strip()
    lower = force_str.lower().replace(" ", "")

    # --- Special signals ---

    # Impulse: "impulse", "delta", "3impulse", "2delta"
    m_imp = re.fullmatch(r"([+-]?\d*\.?\d*)\s*(?:impulse|delta)(?:\(t\))?", lower)
    if m_imp:
        coeff = m_imp.group(1)
        F0 = float(coeff) if coeff and coeff not in ("+", "-", "") else 1.0
        if coeff == "-":
            F0 = -1.0
        return {"type": "impulse", "F0": F0}

    # Step: "step", "u(t)", "u", "3step", "5u(t)"
    m_step = re.fullmatch(r"([+-]?\d*\.?\d*)\s*(?:step|u(?:\(t\))?)(?:\(t\))?", lower)
    if m_step:
        coeff = m_step.group(1)
        F0 = float(coeff) if coeff and coeff not in ("+", "-", "") else 1.0
        if coeff == "-":
            F0 = -1.0
        return {"type": "step", "F0": F0}

    # Ramp: "ramp", "t", "3ramp", "2t"
    m_ramp = re.fullmatch(r"([+-]?\d*\.?\d*)\s*(?:ramp|t)", lower)
    if m_ramp:
        coeff = m_ramp.group(1)
        F0 = float(coeff) if coeff and coeff not in ("+", "-", "") else 1.0
        if coeff == "-":
            F0 = -1.0
        return {"type": "ramp", "F0": F0}

    # Exponential: "exp(-2t)", "e^(-3t)", "5exp(-2t)"
    m_exp = re.fullmatch(
        r"([+-]?\d*\.?\d*)\s*(?:exp|e\^?)\s*\(\s*(-?\d*\.?\d+)\s*\*?\s*t\s*\)",
        lower,
    )
    if m_exp:
        coeff = m_exp.group(1)
        F0 = float(coeff) if coeff and coeff not in ("+", "-", "") else 1.0
        if coeff == "-":
            F0 = -1.0
        a = float(m_exp.group(2))
        return {"type": "exponential", "F0": F0, "a": a}

    # Pure constant: just a number like "3", "0.5"
    m_const = re.fullmatch(r"([+-]?\d+\.?\d*)", lower)
    if m_const:
        return {"type": "step", "F0": float(m_const.group(1))}

    # --- Harmonic signals ---
    _TERM_RE = re.compile(
        r"([+-]?\s*(?:\d+(?:\.\d*)?|\.\d+)?)\s*"  # optional coefficient
        r"(cos|sin)\s*"                             # cos or sin
        r"(?:\(\s*([\d.]+|w|omega|ω)\s*(?:\*?\s*t)?\s*\)"  # (freq) or (w*t)
        r"|([\d.]+|w|omega|ω)?)"                   # freq without parens
        r"\s*(?:\*?\s*t)?",                          # optional trailing t
        re.IGNORECASE,
    )

    matches = list(_TERM_RE.finditer(force_str))
    if not matches:
        raise ValueError(f"Cannot parse forcing function: '{force_str}'")

    components = []
    for m in matches:
        raw_coeff = m.group(1).replace(" ", "") if m.group(1) else ""
        func = m.group(2).lower()
        freq_paren = m.group(3)
        freq_bare = m.group(4)

        if raw_coeff in ("", "+"):
            F0 = 1.0
        elif raw_coeff == "-":
            F0 = -1.0
        else:
            F0 = float(raw_coeff)

        raw_freq = (freq_paren or freq_bare or "").strip().lower()
        if raw_freq in ("w", "omega", "ω", ""):
            omega = "w"  # symbolic — solver will handle
        else:
            omega = float(raw_freq)

        components.append({"F0": abs(F0), "omega": omega, "func": func,
                           "phase_deg": 180.0 if F0 < 0 else 0.0})

    if len(components) == 1:
        c = components[0]
        return {
            "type": "harmonic",
            "F0": c["F0"],
            "omega": c["omega"],
            "func": c["func"],
        }
    else:
        return {
            "type": "general_harmonic",
            "components": components,
        }


def _input_with_default(prompt: str, default: str = "") -> str:
    """Prompt user with a default shown in brackets."""
    if default:
        display = f"{prompt} [{default}]: "
    else:
        display = f"{prompt}: "
    val = input(display).strip()
    return val if val else default


def _input_optional_vector(prompt: str) -> list[float] | None:
    """Prompt for an optional vector; return None if skipped."""
    val = input(f"{prompt} (Enter to skip): ").strip()
    if not val:
        return None
    return parse_vector(val)


# ---------------------------------------------------------------------------
# Display / save helpers
# ---------------------------------------------------------------------------


def display_result(result: SolverResult, mode: str = "exam") -> None:
    """Render and print a SolverResult.

    mode: "exam" (full steps) or "answer" (concise answer with reasoning).
    """
    print()
    if mode == "answer":
        print(ReportEngine.render_exam_answer(result))
    else:
        print(ReportEngine.render_exam_style(result))
    print()


def ask_plot(result: SolverResult) -> None:
    """그래프를 보여줄지 사용자에게 묻고, 원하면 적절한 플롯을 생성한다."""
    try:
        answer = input("그래프를 보시겠습니까? (y/n): ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print()
        return

    if answer not in ("y", "yes", "ㅛ"):
        return

    # Lazy import to avoid loading matplotlib at startup
    try:
        import matplotlib.pyplot as plt
        from ..core.visualizer import Visualizer
    except ImportError as exc:
        print(f"  [오류] 시각화 모듈을 불러올 수 없습니다: {exc}")
        return

    viz = Visualizer()
    # Ask whether to save before plotting
    try:
        save_ans = input("  그래프를 파일로 저장하시겠습니까? (y/n): ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        save_ans = "n"

    do_save = save_ans in ("y", "yes", "ㅛ")
    figs = viz.auto_plot(result, save=do_save)

    if not figs:
        print("  이 결과 유형에 대한 그래프가 없습니다.")
        return

    plt.show()
    print()


def ask_save(result: SolverResult) -> None:
    """Ask the user whether to save the result to a file."""
    try:
        answer = input("결과를 파일로 저장하시겠습니까? (y/n): ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print()
        return

    if answer in ("y", "yes", "ㅛ"):
        _REPORTS_DIR.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        ptype = result.problem_type.replace(" ", "_").replace("/", "_")[:40]
        md_path = _REPORTS_DIR / f"{timestamp}_{ptype}.md"
        txt_path = _REPORTS_DIR / f"{timestamp}_{ptype}.txt"
        try:
            ReportEngine.save_markdown(result, md_path)
            ReportEngine.save_exam_style(result, txt_path)
            print(f"  -> Markdown: {md_path}")
            print(f"  -> Text:     {txt_path}")
        except Exception as e:
            print(f"  저장 실패: {e}")
    print()


def _pause() -> None:
    """Pause before returning to menu."""
    try:
        input("Enter를 눌러 메뉴로 돌아가기...")
    except (EOFError, KeyboardInterrupt):
        pass


# ---------------------------------------------------------------------------
# Menu handlers
# ---------------------------------------------------------------------------


def _handle_concept_db() -> None:
    """1. T/F 개념 판별 (v2: 구조 분석 + 충돌 검사)."""
    from ..core.concept_db import ConceptDBSolver

    print("\n=== T/F 개념 판별 (v2) ===")
    print("T/F 문장을 입력하세요. 솔버가 구조 분석 후 판정합니다.")
    print("예: 'gyroscopic force does no work'")
    print("    'convolution is used for both free and forced vibration'")
    print()

    statement = input("문장: ").strip()
    if not statement:
        print("  입력이 없습니다.")
        return

    extra_kw_str = input("추가 키워드 (쉼표 구분, Enter to skip): ").strip()
    extra_keywords = [k.strip() for k in extra_kw_str.split(",") if k.strip()] if extra_kw_str else []

    topic = input("토픽 필터 (Enter to skip): ").strip() or None

    solver = ConceptDBSolver()
    params: dict[str, Any] = {"statement": statement}
    if extra_keywords:
        params["keywords"] = extra_keywords
    if topic:
        params["topic"] = topic

    result = solver.solve(params)

    # --- v2 advisory display ---
    fa = result.final_answer
    print()
    print("=" * 55)
    if fa.get("advisory"):
        print(fa["advisory"])
    print("-" * 55)
    confidence = fa.get("confidence", "")
    adj_verdict = fa.get("adjusted_verdict", fa.get("verdict", ""))
    conflicts = fa.get("conflicts", [])
    if confidence == "CONFLICT_DETECTED":
        print(f"  ** 추정 판정: {adj_verdict} (충돌 {len(conflicts)}건 감지) **")
    elif confidence == "CAUTION":
        print(f"  ** 추정 판정: {adj_verdict} (주의 필요) **")
    else:
        print(f"  ** 판정: {adj_verdict} (confidence: {confidence}) **")
    print("=" * 55)
    print()

    # Also show full steps via standard display
    display_result(result)
    ask_save(result)


def _handle_frf() -> None:
    """2. FRF 분해."""
    from ..templates.frf_decompose import FRFDecomposeSolver

    print("\n=== FRF 분해 (전달함수) ===")
    print("두 가지 입력 방식 중 선택하세요:")
    print("  a) ODE 문자열  (예: x''' + 3x'' + 6x' + 8x = f(t))")
    print("  b) 계수 리스트 (예: 1,3,6,8)")
    print()

    mode = _input_with_default("입력 방식 (a/b)", "a").lower()

    params: dict[str, Any] = {}

    if mode == "b":
        coeff_str = input("계수 입력 (고차->저차, 쉼표 구분): ").strip()
        if not coeff_str:
            print("  입력이 없습니다.")
            return
        coefficients = [float(c.strip()) for c in coeff_str.split(",")]
        params["coefficients"] = coefficients
    else:
        ode_str = input("ODE 입력: ").strip()
        if not ode_str:
            print("  입력이 없습니다.")
            return
        params["ode_string"] = ode_str

    # Optional numerator
    num_str = input("분자 계수 (Enter for [1]): ").strip()
    if num_str:
        params["numerator"] = [float(c.strip()) for c in num_str.split(",")]

    # Optional forcing function
    print()
    print("강제 입력 f(t) (Enter to skip):")
    print("  조화: cos, 3cos2, 5sin(3), 2cos(1)+3sin(4)")
    print("  심볼릭: 10sin, 10sinw, 5cos(w) → x(t) = .../sqrt(...) 형태")
    print("  특수: impulse, step, 3step, ramp, exp(-2t)")
    force_str = input("f(t) = ").strip()
    if force_str:
        params["forcing"] = _parse_forcing_string(force_str)

    solver = FRFDecomposeSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)
    ask_plot(result)


def _handle_state_space() -> None:
    """3. State-Space / 전이행렬 (상태방정식)."""
    from ..core.state_space import StateSpaceSolver

    print("\n=== State-Space / 전이행렬 (상태방정식) ===")
    print("두 가지 입력 방식 중 선택하세요:")
    print("  a) 상태행렬 직접 입력 (A, B, x0)")
    print("  b) 2차 시스템 변환 (M, C, K → state-space)")
    print()

    mode = _input_with_default("입력 방식 (a/b)", "a").lower()

    params: dict[str, Any] = {}

    if mode == "b":
        # Second-order mode
        params["mode"] = "second_order"
        print("\nMq̈ + Cq̇ + Kq = B_force·f(t)")
        print("행렬 입력: 행은 ';', 열은 ',' 구분 (예: 1,0;0,1)")
        print()

        M_str = input("M (질량행렬): ").strip()
        if not M_str:
            print("  입력이 없습니다.")
            return
        params["M"] = parse_matrix(M_str)

        K_str = input("K (강성행렬): ").strip()
        if not K_str:
            print("  입력이 없습니다.")
            return
        params["K"] = parse_matrix(K_str)

        C_str = input("C (감쇠행렬, Enter for 0): ").strip()
        if C_str:
            params["C"] = parse_matrix(C_str)

        B_str = input("B_force (입력행렬, Enter to skip): ").strip()
        if B_str:
            params["B_force"] = parse_matrix(B_str)

        q0 = _input_optional_vector("초기 변위 q(0) (예: 1,0)")
        qdot0 = _input_optional_vector("초기 속도 qdot(0) (예: 0,0)")
        if q0 is not None:
            params["q0"] = q0
        if qdot0 is not None:
            params["qdot0"] = qdot0
    else:
        # Direct state-matrix mode
        params["mode"] = "state_matrix"
        print("\nẋ = Ax + Bu")
        print("행렬 입력: 행은 ';', 열은 ',' 구분 (예: 0,1;-4,-2)")
        print()

        A_str = input("A (상태행렬): ").strip()
        if not A_str:
            print("  입력이 없습니다.")
            return
        params["A"] = parse_matrix(A_str)

        B_str = input("B (입력행렬, Enter to skip): ").strip()
        if B_str:
            params["B"] = parse_matrix(B_str)

        x0 = _input_optional_vector("초기 상태 x(0) (예: 1,0)")
        if x0 is not None:
            params["x0"] = x0

    # Optional forcing
    print()
    print("강제 입력 f(t) (Enter to skip):")
    print("  조화: cos, 3cos2, 5sin(3), 2cos(1)+3sin(4)")
    print("  특수: impulse, step, 3step, ramp, exp(-2t)")
    force_str = input("f(t) = ").strip()
    if force_str:
        params["forcing"] = _parse_forcing_string(force_str)

    # Optional time span
    t_end_str = input("시뮬레이션 종료 시간 (Enter for 10): ").strip()
    if t_end_str:
        params["t_span"] = [0, float(t_end_str)]

    # Optional evaluation times
    t_eval_str = input("전이행렬 평가 시점 (쉼표 구분, Enter to skip): ").strip()
    if t_eval_str:
        params["t_eval"] = [float(t.strip()) for t in t_eval_str.split(",")]

    solver = StateSpaceSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)
    ask_plot(result)


def _handle_inverse_laplace() -> None:
    """4. Inverse Laplace Transform."""
    from ..core.inverse_laplace import InverseLaplaceSolver

    print("\n=== Inverse Laplace Transform (역 라플라스) ===")
    print("F(s) = N(s)/D(s)의 역 라플라스 변환을 계산합니다.")
    print("두 가지 입력 방식 중 선택하세요:")
    print("  a) 수식 문자열 (예: 1/(s**2 + 4*s + 9))")
    print("  b) 분자/분모 계수 리스트")
    print()

    mode = _input_with_default("입력 방식 (a/b)", "a").lower()

    params: dict[str, Any] = {}

    if mode == "b":
        num_str = input("분자 계수 (고차→저차, 쉼표 구분, 예: 1): ").strip()
        if not num_str:
            print("  입력이 없습니다.")
            return
        params["numerator"] = [float(c.strip()) for c in num_str.split(",")]

        den_str = input("분모 계수 (고차→저차, 쉼표 구분, 예: 1,4,9): ").strip()
        if not den_str:
            print("  입력이 없습니다.")
            return
        params["denominator"] = [float(c.strip()) for c in den_str.split(",")]
    else:
        expr_str = input("F(s) = ").strip()
        if not expr_str:
            print("  입력이 없습니다.")
            return
        params["expression"] = expr_str

    solver = InverseLaplaceSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)


def _handle_lagrange() -> None:
    """5. Lagrange -> equilibrium -> linearization -> stability."""
    from ..core.lagrange import LagrangeSolver

    print("\n=== Lagrange 방정식 분석 ===")
    print("운동에너지 T와 위치에너지 V를 sympy 표현식으로 입력하세요.")
    print("  Rational(1,2)로 분수를, ** 로 거듭제곱을 표현합니다.")
    print()

    dof_str = _input_with_default("자유도 (DOF)", "1")
    n_dof = int(dof_str) if dof_str.isdigit() else 1

    T_expr = input("T (운동에너지) = ").strip()
    if not T_expr:
        print("  입력이 없습니다.")
        return

    V_expr = _input_with_default("V (위치에너지)", "0")

    param_str = input("매개변수 이름 (쉼표 구분, 예: m,g,R,Omega): ").strip()
    parameters = [p.strip() for p in param_str.split(",") if p.strip()] if param_str else []

    param_values = None
    if parameters:
        pv_str = input("매개변수 수치값 (예: m=1,g=9.81 / Enter to skip): ").strip()
        if pv_str:
            param_values = {}
            for item in pv_str.split(","):
                item = item.strip()
                if "=" in item:
                    k, v = item.split("=", 1)
                    param_values[k.strip()] = float(v.strip())

    solver = LagrangeSolver()

    if n_dof > 1:
        # Multi-DOF input
        coords = []
        for i in range(1, n_dof + 1):
            q_name = _input_with_default(f"좌표 {i} 이름", f"q{i}")
            qd_name = _input_with_default(f"속도 {i} 이름", q_name + "dot")
            coords.append((q_name, qd_name))

        Qnc_exprs = {}
        qnc_str = input("비보존력 (예: q1=-c*q1dot,q2=0 / Enter to skip): ").strip()
        if qnc_str:
            for item in qnc_str.split(","):
                item = item.strip()
                if "=" in item:
                    k, v = item.split("=", 1)
                    Qnc_exprs[k.strip()] = v.strip()

        # Multi-DOF: 분석 범위 선택
        print()
        print("분석 범위 선택:")
        print("  a) EOM만 유도 (빠름)")
        print("  b) EOM + 평형점 + 선형화 + 안정성 (느릴 수 있음)")
        scope = _input_with_default("범위 (a/b)", "a").lower()

        if scope == "b":
            find = ["eom", "equilibrium", "linearization", "stability"]
        else:
            find = ["eom"]

        params: dict[str, Any] = {
            "T_expr": T_expr,
            "V_expr": V_expr,
            "coords": coords,
            "parameters": parameters,
            "Qnc_exprs": Qnc_exprs,
            "find": find,
        }
    else:
        # 1-DOF input
        coord = _input_with_default("일반화 좌표 이름", "q")
        coord_dot = _input_with_default("일반화 속도 이름", coord + "dot")
        Qnc_expr = _input_with_default("비보존력 Q_nc (Enter for 0)", "0")

        params = {
            "T_expr": T_expr,
            "V_expr": V_expr,
            "coord": coord,
            "coord_dot": coord_dot,
            "parameters": parameters,
            "Qnc_expr": Qnc_expr,
        }

    if param_values:
        params["parameter_values"] = param_values

    result = solver.solve(params)
    display_result(result)
    ask_save(result)


def _handle_modal() -> None:
    """4. 2-DOF Modal Analysis."""
    from ..core.modal import ModalSolver

    print("\n=== 2-DOF Modal Analysis ===")
    print("M, K 행렬을 입력하세요.")
    print("  형식: 행은 ';' 으로, 열은 ',' 으로 구분")
    print("  예: 9,0;0,9")
    print()

    M_str = input("M (질량행렬): ").strip()
    if not M_str:
        print("  입력이 없습니다.")
        return
    M = parse_matrix(M_str)

    K_str = input("K (강성행렬): ").strip()
    if not K_str:
        print("  입력이 없습니다.")
        return
    K = parse_matrix(K_str)

    q0 = _input_optional_vector("초기 변위 q(0) (예: 1,0)")
    qdot0 = _input_optional_vector("초기 속도 qdot(0) (예: 0,0)")

    params: dict[str, Any] = {"M": M, "K": K}
    if q0 is not None:
        params["initial_q"] = q0
    if qdot0 is not None:
        params["initial_qdot"] = qdot0

    solver = ModalSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)
    ask_plot(result)


def _handle_damping() -> None:
    """5. Proportional Damping Response."""
    from ..core.damping import DampingSolver

    print("\n=== Proportional Damping Response ===")
    print("M, K 행렬과 감쇠 정보를 입력하세요.")
    print()

    M_str = input("M (질량행렬, 예: 9,0;0,9): ").strip()
    if not M_str:
        print("  입력이 없습니다.")
        return
    M = parse_matrix(M_str)

    K_str = input("K (강성행렬, 예: 72,-36;-36,72): ").strip()
    if not K_str:
        print("  입력이 없습니다.")
        return
    K = parse_matrix(K_str)

    print()
    print("감쇠 입력 방식:")
    print("  a) C 행렬 직접 입력")
    print("  b) alpha, beta (Rayleigh damping: C = alpha*M + beta*K)")
    print("  c) 모달 감쇠비 zeta 직접 입력")
    damp_mode = _input_with_default("감쇠 입력 방식 (a/b/c)", "b").lower()

    params: dict[str, Any] = {"M": M, "K": K}

    if damp_mode == "a":
        C_str = input("C (감쇠행렬): ").strip()
        if not C_str:
            print("  입력이 없습니다.")
            return
        params["C"] = parse_matrix(C_str)
    elif damp_mode == "c":
        zeta_str = input("모달 감쇠비 (쉼표 구분, 예: 0.05,0.1): ").strip()
        if not zeta_str:
            print("  입력이 없습니다.")
            return
        params["zeta"] = parse_vector(zeta_str)
    else:
        alpha_str = input("alpha: ").strip()
        beta_str = input("beta: ").strip()
        if not alpha_str or not beta_str:
            print("  alpha와 beta를 모두 입력하세요.")
            return
        params["alpha"] = float(alpha_str)
        params["beta"] = float(beta_str)

    q0 = _input_optional_vector("초기 변위 q(0)")
    qdot0 = _input_optional_vector("초기 속도 qdot(0)")
    if q0 is not None:
        params["initial_q"] = q0
    if qdot0 is not None:
        params["initial_qdot"] = qdot0

    solver = DampingSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)
    ask_plot(result)


def _handle_auto_route() -> None:
    """a. 자동 분류 — 문제 텍스트 입력 -> 유형 판별 -> 해당 메뉴 호출."""
    from ..core.router import classify, get_menu_number

    print("\n=== 자동 문제 분류 ===")
    print("문제 텍스트 또는 키워드를 입력하세요.")
    print("예: 'x\\'\\'\\' + 3x\\'\\'  + 6x\\' + 8x = f(t)'")
    print("    'T = 1/2*m*rdot^2, V = -mgR*cos(theta)'")
    print("    '이 문장이 참인지 거짓인지 판별'")
    print()

    text = input("입력: ").strip()
    if not text:
        print("  입력이 없습니다.")
        return

    result = classify(text)

    print(f"\n  분류 결과: [{result.route}] {result.label}")
    print(f"  신뢰도: {result.confidence:.0%}")
    if result.matched_keywords:
        kws = [k for k in result.matched_keywords if not k.startswith("pattern:")]
        if kws:
            print(f"  매칭 키워드: {', '.join(kws[:5])}")
    if result.suggestions:
        print(f"  대안: {', '.join(result.suggestions)}")
    print()

    if result.confidence < 0.1:
        print("  분류 신뢰도가 너무 낮습니다. 메뉴에서 직접 선택해주세요.")
        return

    menu_num = get_menu_number(result.route)
    if menu_num:
        try:
            answer = input(f"  [{result.label}] 메뉴를 실행하시겠습니까? (y/n) [y]: ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print()
            return
        if answer in ("", "y", "yes", "ㅛ"):
            handler_entry = _HANDLERS.get(menu_num)
            if handler_entry:
                _, handler_fn = handler_entry
                if handler_fn:
                    handler_fn()
    else:
        print("  해당 라우트에 연결된 메뉴가 없습니다.")


def _handle_extended() -> None:
    """6. 확장 (Gyroscopic / Nonconservative)."""
    from ..core.extended import ExtendedSolver

    print("\n=== 확장 시스템 (Gyroscopic / Nonconservative) ===")
    print("M q'' + (C+G) q' + (K+H) q = 0")
    print("  M: 질량행렬 (대칭, 양정치)")
    print("  K: 강성행렬 (대칭)")
    print("  G: 자이로 행렬 (반대칭, 선택)")
    print("  H: 비보존력 행렬 (반대칭, 선택)")
    print("  C: 감쇠행렬 (대칭, 선택)")
    print("행렬 입력: 행은 ';', 열은 ',' 구분 (예: 1,0;0,1)")
    print()

    M_str = input("M (질량행렬): ").strip()
    if not M_str:
        print("  입력이 없습니다.")
        return
    M = parse_matrix(M_str)

    K_str = input("K (강성행렬): ").strip()
    if not K_str:
        print("  입력이 없습니다.")
        return
    K = parse_matrix(K_str)

    G_str = input("G (자이로 행렬, Enter to skip): ").strip()
    G = parse_matrix(G_str) if G_str else None

    H_str = input("H (비보존력 행렬, Enter to skip): ").strip()
    H = parse_matrix(H_str) if H_str else None

    C_str = input("C (감쇠행렬, Enter to skip): ").strip()
    C = parse_matrix(C_str) if C_str else None

    params: dict[str, Any] = {"M": M, "K": K}
    if G:
        params["G"] = G
    if H:
        params["H"] = H
    if C:
        params["C"] = C

    # Optional initial conditions
    q0 = _input_optional_vector("초기 변위 q(0) (Enter to skip)")
    qdot0 = _input_optional_vector("초기 속도 qdot(0) (Enter to skip)")
    if q0:
        params["initial_q"] = q0
    if qdot0:
        params["initial_qdot"] = qdot0

    solver = ExtendedSolver()
    result = solver.solve(params)
    display_result(result)

    # Offer plot
    ask_plot(result)
    ask_save(result)


def _handle_base_excitation() -> None:
    """9. Base Excitation (지반 가진)."""
    from ..core.base_excitation import BaseExcitationSolver, parse_base_excitation_eq

    print("\n=== Base Excitation (지반 가진) ===")
    print("입력 방식:")
    print("  a) 식 입력  (예: 2x''+3x'+5x = 3y'+5y, y=0.05sin8t)")
    print("  b) 파라미터 개별 입력 (M, c, k, A, w)")
    print("  c) MDOF (다자유도 행렬)")
    print()

    mode = _input_with_default("입력 방식 (a/b/c)", "a").lower()
    params: dict[str, Any] = {}

    if mode == "c":
        # --- MDOF ---
        print("\nMx'' + Cx' + Kx = Ci y' + Ki y")
        print("행렬: 행은 ';', 열은 ',' 구분")
        print()

        M_str = input("M: ").strip()
        if not M_str:
            print("  입력이 없습니다.")
            return
        params["M"] = parse_matrix(M_str)

        K_str = input("K: ").strip()
        if not K_str:
            print("  입력이 없습니다.")
            return
        params["K"] = parse_matrix(K_str)

        params["n_dof"] = len(params["M"])

        C_str = input("C (Enter for 0): ").strip()
        if C_str:
            params["C"] = parse_matrix(C_str)

        iota_str = input("influence vector (Enter for all 1s): ").strip()
        if iota_str:
            params["influence_vector"] = parse_vector(iota_str)

        params["A"] = float(_input_with_default("A", "1"))
        params["omega"] = float(_input_with_default("w (rad/s)", "1"))

    elif mode == "a":
        # --- 식 입력 ---
        print()
        print("EOM 예시:")
        print("  2x''+3x'+5x = 3y'+5y")
        print("  Mx''+c(x'-y')+k(x-y) = 0")
        print("  x''+2x'+4x = 2y'+4y")
        print()
        eom_str = input("EOM: ").strip()
        if not eom_str:
            print("  입력이 없습니다.")
            return

        print()
        print("가진 예시:  y=0.05sin8t  /  y=Asinwt  /  y(t)=3sin(10t)")
        print()
        exc_str = input("y(t) = ").strip()
        if not exc_str:
            exc_str = "y=Asinwt"
        if not exc_str.lower().startswith("y"):
            exc_str = "y=" + exc_str

        parsed = parse_base_excitation_eq(eom_str, exc_str)
        params["n_dof"] = 1

        print()
        print("  [파싱 결과]")
        labels = {"m": "M", "c": "c", "k": "k", "A": "A", "omega": "w"}
        for key in ["m", "c", "k", "A", "omega"]:
            val = parsed.get(key)
            if val is not None:
                params[key] = val
                print(f"    {labels[key]} = {val}")
            else:
                print(f"    {labels[key]} = (symbolic)")
        print()

    else:
        # --- 파라미터 개별 입력 ---
        params["n_dof"] = 1
        print("\n수치 입력시 유도 + 평가, Enter시 symbolic 유도만.")
        print()

        for key, label in [("m", "M"), ("k", "k"), ("c", "c"), ("A", "A"), ("omega", "w")]:
            val_str = input(f"{label} (Enter=symbolic): ").strip()
            if val_str:
                params[key] = float(val_str)

    solver = BaseExcitationSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)


def _handle_fourier() -> None:
    """10. Fourier 급수 (주기 가진 응답)."""
    from ..core.fourier import FourierSolver

    print("\n=== Fourier 급수 (주기 가진 응답) ===")
    print("주기 함수를 Fourier 급수로 전개하고, 각 조화 성분의")
    print("정상상태 응답을 중첩하여 전체 응답을 구합니다.")
    print()

    # 시스템 차수
    order_str = _input_with_default("시스템 차수 (1: cx'+kx=f, 2: mx''+cx'+kx=f)", "1")
    system_order = int(order_str)

    # 파형 선택
    print("\n파형 종류:")
    print("  1) square_wave  (구형파)")
    print("  2) sawtooth     (톱니파)")
    print("  3) triangle     (삼각파)")
    print("  4) custom       (사용자 정의 계수)")
    wave_choice = _input_with_default("파형 (1-4)", "1")
    waveform_map = {"1": "square_wave", "2": "sawtooth", "3": "triangle", "4": "custom"}
    waveform = waveform_map.get(wave_choice, "square_wave")

    params: dict[str, Any] = {
        "system_order": system_order,
        "waveform": waveform,
    }

    # 시스템 파라미터
    print("\n수치 입력시 유도 + 평가, Enter시 symbolic 유도만.")
    print()

    if system_order >= 2:
        m_str = input("m (Enter=symbolic): ").strip()
        params["m"] = m_str if m_str else "m"

    c_str = input("c (Enter=symbolic): ").strip()
    params["c"] = float(c_str) if c_str else "c"

    k_str = input("k (Enter=symbolic): ").strip()
    params["k"] = float(k_str) if k_str else "k"

    # 가진 파라미터
    f0_str = input("f0 (가진력 진폭, Enter=symbolic): ").strip()
    params["f0"] = float(f0_str) if f0_str else "f0"

    T_str = input("T (주기, Enter=symbolic): ").strip()
    params["T"] = float(T_str) if T_str else "T"

    N_str = _input_with_default("Fourier 항 수 N", "10")
    params["num_terms"] = int(N_str)

    # custom 파형일 때 계수 입력
    if waveform == "custom":
        a0_str = input("a0: ").strip()
        params["a0"] = float(a0_str) if a0_str else 0
        an_str = input("an (쉼표 구분, 예: 1,0,0.5): ").strip()
        if an_str:
            params["an"] = [float(x) for x in an_str.split(",")]
        bn_str = input("bn (쉼표 구분): ").strip()
        if bn_str:
            params["bn"] = [float(x) for x in bn_str.split(",")]

    solver = FourierSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)


def _handle_convolution() -> None:
    """11. Convolution Integral (Duhamel 적분)."""
    from ..core.convolution import ConvolutionSolver

    print("\n=== Convolution Integral (Duhamel 적분) ===")
    print("x(t) = ∫₀ᵗ h(t-τ)·f(τ) dτ  (zero initial conditions)")
    print()

    # 시스템 차수
    order_str = _input_with_default("시스템 차수 (1: cx'+kx=f, 2: mx''+cx'+kx=f)", "2")
    system_order = int(order_str)

    params: dict[str, Any] = {"system_order": system_order}

    # 시스템 파라미터
    print("\n수치 입력시 유도 + 평가, Enter시 symbolic 유도만.")
    print()

    if system_order >= 2:
        m_str = input("m (Enter=symbolic): ").strip()
        params["m"] = float(m_str) if m_str else "m"

    c_str = input("c (Enter=symbolic): ").strip()
    params["c"] = float(c_str) if c_str else "c"

    k_str = input("k (Enter=symbolic): ").strip()
    params["k"] = float(k_str) if k_str else "k"

    # 가진 함수 입력
    print("\n가진 함수 f(t) 입력 방식:")
    print("  a) 프리셋 (step, ramp, exp, exp_truncated, harmonic, pulse)")
    print("  b) 구간별 직접 입력 (piecewise)")
    print()

    mode = _input_with_default("입력 방식 (a/b)", "a").lower()

    if mode == "b":
        # --- Piecewise 입력 ---
        params["input_mode"] = "piecewise"
        print("\n구간별 함수를 입력하세요. (빈 줄 입력시 종료)")
        print("예시:")
        print("  구간 1: expr=F0*exp(-a*t), start=0, end=t1")
        print("  구간 2: expr=0, start=t1, end=inf")
        print()
        print("사용 가능 변수: t, F0, f0, a, b, t1, t2, T, omega")
        print()

        pieces = []
        idx = 1
        while True:
            print(f"--- 구간 {idx} ---")
            expr_str = input("  f(t) = ").strip()
            if not expr_str:
                break
            start_str = input("  시작 (default=0): ").strip() or "0"
            end_str = input("  끝 (default=inf): ").strip() or "inf"
            pieces.append({
                "expr": expr_str,
                "start": start_str,
                "end": end_str,
            })
            idx += 1

        if not pieces:
            print("  입력이 없습니다.")
            return
        params["pieces"] = pieces

    else:
        # --- 프리셋 ---
        params["input_mode"] = "preset"
        print("\n프리셋 종류:")
        print("  1) step             F0·u(t)")
        print("  2) ramp             F0·t")
        print("  3) exponential      F0·e^(-at)")
        print("  4) exp_truncated    F0·e^(-at), 0<t<t1 / 0, t>t1")
        print("  5) harmonic         F0·sin(ωt) or cos(ωt)")
        print("  6) rectangular_pulse  F0, 0<t<t1 / 0, t>t1")
        print()
        preset_map = {
            "1": "step", "2": "ramp", "3": "exponential",
            "4": "exp_truncated", "5": "harmonic", "6": "rectangular_pulse",
        }
        preset_choice = _input_with_default("선택 (1-6)", "4")
        force_type = preset_map.get(preset_choice, "step")
        params["force_type"] = force_type

        F0_str = input("F0 (Enter=symbolic): ").strip()
        params["F0"] = float(F0_str) if F0_str else "F0"

        if force_type in ("exponential", "exp_truncated"):
            a_str = input("a (감쇠율, Enter=symbolic): ").strip()
            params["a"] = float(a_str) if a_str else "a"

        if force_type in ("exp_truncated", "rectangular_pulse"):
            t1_str = input("t1 (구간 끝, Enter=symbolic): ").strip()
            params["t1"] = float(t1_str) if t1_str else "t1"

        if force_type == "harmonic":
            omega_str = input("omega (Enter=symbolic): ").strip()
            params["omega"] = float(omega_str) if omega_str else "omega"
            func_str = _input_with_default("sin or cos", "sin")
            params["func"] = func_str

    solver = ConvolutionSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)


def _handle_rotating_hoop() -> None:
    """7. Rotating Hoop template."""
    from ..templates.rotating_hoop import RotatingHoopSolver

    print("\n=== 빈출 템플릿: Rotating Hoop ===")
    print("회전하는 후프 위의 비드 문제입니다.")
    print("숫자를 입력하거나, 'symbolic' / Enter로 기호 해석합니다.")
    print()

    R = parse_float_or_symbol(_input_with_default("R (후프 반지름, m)", "symbolic"))
    m = parse_float_or_symbol(_input_with_default("m (질점 질량, kg)", "symbolic"))
    Omega = parse_float_or_symbol(_input_with_default("Omega (각속도, rad/s)", "symbolic"))
    g = parse_float_or_symbol(_input_with_default("g (중력가속도, m/s^2)", "9.81"))

    params: dict[str, Any] = {"R": R, "m": m, "Omega": Omega, "g": g}

    solver = RotatingHoopSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)


def _handle_rotating_triangle() -> None:
    """8. Rotating Triangle template."""
    from ..templates.rotating_triangle import RotatingTriangleSolver

    print("\n=== 빈출 템플릿: Rotating Triangle ===")
    print("회전 삼각형 프레임 위의 질점 문제입니다.")
    print()

    orient = _input_with_default("orientation (midterm_2025 / simple)", "midterm_2025")

    m = parse_float_or_symbol(_input_with_default("m (질점 질량, kg)", "symbolic"))
    Omega = parse_float_or_symbol(_input_with_default("Omega (각속도, rad/s)", "symbolic"))
    g = parse_float_or_symbol(_input_with_default("g (중력가속도, m/s^2)", "9.81"))

    params: dict[str, Any] = {"m": m, "Omega": Omega, "g": g, "orientation": orient}

    if orient == "simple":
        angle = _input_with_default("rod angle from horizontal (deg)", "30")
        params["orientation_angle"] = float(angle)

    solver = RotatingTriangleSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)


def _handle_two_dof_chain() -> None:
    """9. 2-DOF Chain System template."""
    from ..templates.two_dof_chain import TwoDOFChainSolver

    print("\n=== 빈출 템플릿: 2-DOF Chain System ===")
    print("  wall -- k1 -- [m1] -- k2 -- [m2] -- k3 -- wall")
    print()

    m_str = input("질량 [m1, m2] (쉼표 구분, 예: 9,9): ").strip()
    if not m_str:
        print("  입력이 없습니다.")
        return
    masses = [float(x.strip()) for x in m_str.split(",")]

    k_str = input("스프링 [k1, k2, k3] (쉼표 구분, 예: 36,36,36 / k3=0이면 36,36): ").strip()
    if not k_str:
        print("  입력이 없습니다.")
        return
    springs = [float(x.strip()) for x in k_str.split(",")]

    topology = _input_with_default("topology (chain / free_free)", "chain")

    params: dict[str, Any] = {
        "masses": masses,
        "springs": springs,
        "topology": topology,
    }

    # Optional initial conditions
    ic_str = input("초기조건 입력? (y/n) [n]: ").strip().lower()
    if ic_str in ("y", "yes"):
        x0 = _input_optional_vector("초기 변위 x(0) (예: 1,0)")
        v0 = _input_optional_vector("초기 속도 v(0) (예: 0,0)")
        ic: dict[str, Any] = {}
        if x0 is not None:
            ic["x0"] = x0
        if v0 is not None:
            ic["v0"] = v0
        if ic:
            params["initial_conditions"] = ic

    zeta_str = input("모달 감쇠비 zeta (Enter to skip): ").strip()
    if zeta_str:
        params["damping_ratio"] = float(zeta_str)

    solver = TwoDOFChainSolver()
    result = solver.solve(params)
    display_result(result)
    ask_save(result)


# ---------------------------------------------------------------------------
# Dispatch table
# ---------------------------------------------------------------------------

_HANDLERS: dict[str, tuple[str, Any]] = {
    "a": ("자동 분류", _handle_auto_route),
    "A": ("자동 분류", _handle_auto_route),
    "1": ("T/F 개념 판별", _handle_concept_db),
    "2": ("FRF 분해 (전달함수)", _handle_frf),
    "3": ("State-Space / 전이행렬", _handle_state_space),
    "4": ("Inverse Laplace Transform", _handle_inverse_laplace),
    "5": ("Lagrange 분석", _handle_lagrange),
    "6": ("N-DOF Modal Analysis", _handle_modal),
    "7": ("Proportional Damping", _handle_damping),
    "8": ("확장 (Gyroscopic)", _handle_extended),
    "9": ("Base Excitation", _handle_base_excitation),
    "10": ("Fourier 급수", _handle_fourier),
    "11": ("Convolution Integral", _handle_convolution),
    "12": ("Rotating Hoop", _handle_rotating_hoop),
    "13": ("Rotating Triangle", _handle_rotating_triangle),
    "14": ("2-DOF Chain System", _handle_two_dof_chain),
    "h": ("도움말", None),
    "H": ("도움말", None),
}


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------


def handle_choice(choice: str) -> bool:
    """Process a single menu choice. Returns False to exit."""
    if choice == "0":
        print("\n종료합니다. 시험 화이팅!")
        return False

    if choice.lower() == "h":
        print(HELP_TEXT)
        return True

    handler_entry = _HANDLERS.get(choice)
    if handler_entry is None:
        print(f"\n  잘못된 입력: '{choice}'. a, 0~14, 또는 h 를 입력하세요.\n")
        return True

    _, handler_fn = handler_entry
    if handler_fn is None:
        print(HELP_TEXT)
        return True

    try:
        handler_fn()
    except KeyboardInterrupt:
        print("\n  (취소됨)")
    except Exception as e:
        print(f"\n  [오류] {type(e).__name__}: {e}")
        print("  상세 정보가 필요하면 아래 traceback을 확인하세요:")
        traceback.print_exc()
        print()

    return True


def main() -> None:
    """CLI main loop."""
    print(BANNER)

    while True:
        print(MENU)
        try:
            choice = input("선택 >> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n\n종료합니다. 시험 화이팅!")
            break

        if not choice:
            continue

        if not handle_choice(choice):
            break


if __name__ == "__main__":
    main()
