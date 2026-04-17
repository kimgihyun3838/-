# ME551 선형진동공학 시험용 Solver 개발 계획서 (plan_ver1)

## 1. 목적

이 문서는 **ME551 선형진동공학 중간고사 전용 문제 풀이 프로그램**을 만들기 위한 1차 개발 계획서이다.

핵심 목표는 다음과 같다.

1. **시험 빈출 문제를 빠르게 분류**한다.
2. 문제 유형별로 **고정된 공식/정리/알고리즘**을 적용한다.
3. 숫자 결과만이 아니라 **서술형 풀이 단계**까지 자동 생성한다.
4. **인터넷 없이 오프라인 단독 실행**이 가능하도록 만든다.

즉, 이 프로그램의 정체는 일반적인 수학 툴이나 PDF 챗봇이 아니라,

> **ME551 중간고사 전용 solver + 풀이 리포트 생성기**

여야 한다.

---

## 2. 전제와 범위

강의 소개 자료 기준으로 ME551은 **이론 + 계산(modeling & solutions based on linear algebra)** 중심 과목이며, 중간고사는 사실상 **Chap. 1~4** 범위를 중심으로 구성된다. 또한 최근 중간고사 문제들은 다음 유형이 반복된다.

- T/F + 이유 설명
- 3차 시스템 FRF 분해
- Lagrange 방정식 → 평형점 → 선형화 → 안정성
- 2-DOF 대칭 고유치 문제 및 modal analysis
- 확장 영역으로 gyroscopic / nonconservative / right-left eigenvector / bi-orthogonality

따라서 v1에서는 범위를 아래와 같이 한정한다.

### v1 핵심 범위
- 선형시스템 기초
- SDOF FRF / impulse / convolution
- Lagrange equation
- equilibrium / perturbation / linearization / stability
- undamped 2-DOF symmetric EVP
- proportional damping modal analysis

### v1.5 확장 범위
- state-space 기반 안정성 판별
- transition matrix
- generalized state equation

### v2 확장 범위
- nonproportional damping
- unsymmetric EVP
- left/right eigenvector
- bi-orthogonality
- gyroscopic complex modes

---

## 3. 개발 철학

### 3.1 Solver-first, Retrieval-second
강의노트 검색 기능은 **개념 확인**용으로만 사용하고, 실제 정답 생성은 반드시 **수식 엔진**이 담당해야 한다.

즉,
- 개념은 local note/reference DB에서 찾고,
- 계산은 symbolic/numeric engine으로 하고,
- 최종 출력은 풀이 템플릿으로 정리한다.

### 3.2 NLP-first로 시작하지 않는다
시험장에서 문제 문장을 자동 파싱하는 기능은 좋지만 안정성이 떨어진다. v1은 다음 구조가 가장 안전하다.

- 사용자가 **문제 유형을 직접 선택**
- 필요한 **파라미터를 직접 입력**
- solver가 결과와 풀이를 생성

즉,

> 자유문장 OCR/해석보다 **구조화 입력 + 전용 solver**가 우선이다.

### 3.3 오프라인 단독 실행이 기본
시험 환경에서 인터넷이 금지되어 있으므로 아래는 피한다.

- OpenAI API
- 웹 검색 의존성
- 클라우드 notebook
- 외부 서버 인증이 필요한 소프트웨어

프로그램은 반드시 **로컬 파일 + 로컬 계산 라이브러리**만으로 동작해야 한다.

---

## 4. 기술 선택

## 권장 선택: Python 우선, MATLAB 보조

### Python을 우선 추천하는 이유
- **SymPy**: symbolic derivation, Lagrange 전개, 평형점 계산, 선형화
- **NumPy / SciPy**: 고유치, 행렬 계산, state-space, 수치응답
- **Matplotlib**: FRF, pole, mode shape 시각화
- **PyYAML / JSON**: 문제 규칙과 개념 DB 관리
- 오프라인 배포가 쉽고 CLI/GUI 둘 다 가능

### MATLAB을 보조 추천하는 이유
- `eig`, `expm`, generalized EVP 검산이 매우 편리함
- 수업 표현과 유사한 수식 실험이 쉬움
- 다만 식 전개 자동화와 템플릿화는 Python 쪽이 더 유연함

### 최종 추천
- **주 개발 언어: Python**
- **검산/레퍼런스: MATLAB**

---

## 5. 프로그램의 최종 형태

프로그램은 아래 3개 층으로 구성하는 것이 좋다.

### 5.1 문제 선택 층
사용자가 아래 중 하나를 선택한다.

1. T/F 개념형
2. FRF 분해형
3. Lagrange-평형-선형화-안정성형
4. 2-DOF modal analysis형
5. proportional damping modal response형
6. 확장: gyroscopic / nonconservative / left-right eigenvector형

### 5.2 계산 엔진 층
선택된 문제 유형에 따라 해당 solver가 계산한다.

### 5.3 풀이 리포트 층
결과를 사람이 제출 가능한 형태로 정리한다.

- Given
- Assumptions
- Governing equations
- Derivation
- Final answer
- Sanity check

---

## 6. 핵심 모듈 구조

## 6.1 concept_db
역할: T/F 문제와 서술형 이유 설명을 담당하는 개념 데이터베이스.

### 담아야 할 내용 예시
- 선형성 = homogeneity + additivity
- arbitrary excitation response:
  - time domain에서는 convolution
  - frequency domain에서는 FRF와 입력의 곱
- p차 MIMO 시스템의 1차 state equation 차원은 `p*n`
- Principle of Virtual Work에서 internal force / ideal constraint force 제거
- structural damping은 harmonic excitation / harmonic response에서만 유효
- real symmetric positive definite matrix의 eigenvalue/eigenvector 성질
- gyroscopic / nonconservative system에서 bi-orthogonality 가능

### 출력 예시
- Verdict: True / False
- Reason: 2~5문장
- Related concepts
- Related lecture tags

---

## 6.2 frf_engine
역할: 1차/2차/3차 선형 시스템 FRF 문제 해결.

### 주요 기능
- ODE를 Laplace domain으로 변환
- 전달함수 `G(s)` 생성
- `G(jω)` 계산
- partial fraction decomposition
- 3차 시스템을 1차 + 2차 subsystem으로 분해
- 1차의 `ω_b`, 2차의 `ω_n`, `ζ` 추출
- characteristic point 계산
- 대략적인 magnitude sketch를 위한 기준점 산출

### 지원 문제 예시
- `x''' + 3x'' + 6x' + 8x = f(t)` 류
- `x''' + 1.2x'' + 4.2x' + 4x = f(t)` 류

### v1 출력 포맷
- transfer function
- 분해 결과
- subsystem parameters
- characteristic points
- sketch guidance

---

## 6.3 lagrange_engine
역할: 1-DOF 또는 저자유도 회전/구속계 문제를 symbolic하게 처리.

### 주요 기능
- generalized coordinate 정의
- 위치/속도 표현식 생성
- kinetic energy `T`
- potential energy `V`
- generalized nonconservative work 또는 force `Q_nc`
- Lagrange equation 생성
- equilibrium point 계산
- perturbation variable 도입
- linearization
- stability classification

### 시험형 템플릿 우선 지원
- rotating hoop
- rotating triangle
- bead-on-constraint
- sliding mass on rotating frame

### 구현 방식
일반 symbolic solver + 시험 빈출 템플릿을 같이 둔다.

- 일반 solver: 식을 유연하게 처리
- 템플릿 solver: 입력을 최소화하고 신뢰도를 높임

---

## 6.4 stability_engine
역할: 평형점 근방의 안정성 판단.

### 지원 방식
1. **선형화 후 eigenvalue 판별**
2. 필요 시 **Liapunov 함수 기반 서술형 설명 보조**

### 출력 예시
- Stable
- Asymptotically stable
- Unstable
- 판정 근거:
  - stiffness sign
  - eigenvalue real part sign
  - damping 유무

---

## 6.5 modal_engine
역할: undamped MDOF와 symmetric eigenvalue problem 처리.

### 주요 기능
- `Mq¨ + Kq = 0` 구성
- generalized EVP `Ku = λMu` 풀이
- 또는 `A v = λ v` 형태로 변환
- eigenvalue/eigenvector 계산
- mass normalization
- orthogonality 검증
- modal matrix `U` 생성
- modal coordinate `η(t)` 구성
- physical response `q(t) = Uη(t)` 복원

### 필수 출력
- `M`, `K`
- eigenvalues
- mode shapes
- `U^T M U = I` 확인
- `U^T K U = diag(λ_i)` 확인
- modal initial condition
- final physical solution

---

## 6.6 damping_engine
역할: 감쇠가 있는 modal response 처리.

### v1 지원 범위
- **proportional damping only**
- `C = αM + βK`인지 검사
- modal decoupling 수행
- 각 모드의 `ζ_r`, `ω_n,r`, `ω_d,r` 계산
- SDOF damped free vibration formula 적용

### v2 지원 범위
- nonproportional damping
- full state-space eigenanalysis
- complex modes

---

## 6.7 state_space_engine
역할: 1차 상태방정식 및 transition matrix 처리.

### 주요 기능
- `x_dot = A x + B f`
- `exp(A t)` 또는 `scipy.linalg.expm`
- homogeneous / forced response
- stability by eigenvalues
- unsymmetric matrix handling

### 활용 시점
- 1차 시스템
- 2차 시스템의 state augmentation
- nonproportional / nonconservative 확장

---

## 6.8 report_engine
역할: 계산 결과를 시험 제출형 풀이로 자동 정리.

### 출력 템플릿
1. Problem type
2. Given
3. Definitions / coordinates
4. Governing equations
5. Derivation steps
6. Simplification / linearization
7. Final answer
8. Interpretation

### 중요 포인트
이 과목은 **이유 설명**이 점수에 직접 연결되므로, 숫자 결과만 보여주면 안 된다.
반드시 **중간식과 논리 흐름**이 같이 나와야 한다.

---

## 7. 입력 방식 설계

## v1 입력 원칙: 구조화 입력

### 예시 1: FRF 문제
```yaml
problem_type: frf_decomposition
ode: "x''' + 3x'' + 6x' + 8x = f(t)"
want:
  - G1
  - G2
  - wb
  - wn
  - zeta
```

### 예시 2: 2-DOF modal problem
```yaml
problem_type: modal_2dof
M: [[9,0],[0,9]]
K: [[72,-36],[-36,72]]
initial_q: [1,1]
initial_qdot: [0,0]
damping:
  type: proportional
  C: "0.2*M + 0.1*K"
```

### 예시 3: rotating hoop / triangle 문제
```yaml
problem_type: rotating_triangle
parameters:
  m: symbolic
  L: symbolic
  g: symbolic
  Omega: symbolic
coordinate: r
asks:
  - eom
  - equilibrium
  - linearization
  - stability
```

---

## 8. UI 설계

### 8.1 1순위: CLI
가장 빠르고 안정적이다.

예시:
```bash
python cli.py
```

실행 후:
- 문제 유형 선택
- 값 입력
- 풀이 출력
- 필요 시 markdown 저장

### 8.2 2순위: 간단한 GUI
시간이 남으면 `tkinter` 또는 `gradio (오프라인 가능 시)` 형태를 고려할 수 있으나,
시험용 안정성만 보면 CLI가 더 낫다.

### 8.3 권장 배포 방식
- Python script 그대로 사용
- 또는 `PyInstaller`로 단일 실행 파일 생성

> 시험 환경에서는 설치 이슈가 적은 형태가 중요하므로,
> **최종 목표는 단일 폴더 혹은 단일 실행 파일**이다.

---

## 9. 문제 유형별 우선순위

### Phase 1 — 문제 taxonomy 정리
문제 유형을 먼저 고정한다.

1. T/F 개념형
2. FRF/전달함수형
3. Lagrange-평형-선형화-안정성형
4. 2-DOF modal analysis형
5. proportional damping response형
6. 확장: nonproportional / gyroscopic / bi-orthogonal형

### Phase 2 — T/F 엔진 완성
가장 빨리 만들 수 있고, 점수 기여도가 높다.

우선 담아야 할 함정:
- convolution vs multiplication
- free vibration vs forced vibration
- state dimension
- structural damping validity
- orthogonality / bi-orthogonality

### Phase 3 — FRF 엔진 완성
최근 기출에서 거의 정형화된 문제다.

필수 기능:
- partial fraction
- first-order + second-order decomposition
- `ω_b`, `ω_n`, `ζ`
- characteristic points

### Phase 4 — Lagrange 템플릿 엔진
회전 hoop / triangle / bead 문제를 템플릿화한다.

필수 기능:
- `T, V` 생성
- EOM
- equilibrium
- linearization
- stability

### Phase 5 — 2-DOF modal engine
필수 기능:
- generalized EVP
- mass normalization
- modal response
- physical reconstruction

### Phase 6 — proportional damping engine
필수 기능:
- `C = αM + βK` 확인
- modal decoupling
- damped SDOF response

### Phase 7 — 확장 엔진
남는 시간에 아래를 확장한다.

- state-space general solver
- gyroscopic systems
- nonconservative systems
- right/left eigenvectors
- bi-orthogonality

---

## 10. 시험장에서의 실제 사용 시나리오

가장 현실적인 사용 흐름은 아래와 같다.

1. 문제를 읽고 유형을 선택한다.
2. 필요한 파라미터를 직접 입력한다.
3. solver가 symbolic / numeric 계산을 수행한다.
4. 프로그램이 최종 답 + 풀이를 출력한다.
5. 사용자는 출력 내용을 바탕으로 답안을 작성한다.

즉,

> “문제 사진 자동해석”보다
> **“내가 구조와 숫자를 넣으면 풀이가 바로 나오는 시스템”**
> 이 훨씬 안전하다.

---

## 11. 테스트 전략

프로그램은 반드시 **기출 회귀 테스트 기반**으로 개발해야 한다.

### 11.1 회귀 테스트 세트 구성
- 2025 midterm 각 문항
- 2024 midterm 각 문항
- 2020~2023 selected problem 중 반복 유형

### 11.2 테스트 기준
- 최종 수치 일치
- 부호 일치
- normalization 일치
- 안정성 판정 일치
- 풀이 단계의 논리 순서 유지

### 11.3 추천 방식
각 기출 문제를 다음처럼 저장한다.

```python
case = {
    "name": "midterm_2025_q2",
    "type": "frf_decomposition",
    "input": {...},
    "expected": {...}
}
```

이후 `pytest`로 자동 검증한다.

---

## 12. 추천 디렉터리 구조

```text
me551_solver/
  core/
    concept_db.py
    frf.py
    lagrange.py
    stability.py
    modal.py
    damping.py
    state_space.py
    report.py
  templates/
    rotating_hoop.py
    rotating_triangle.py
    two_dof_chain.py
    frf_decompose.py
  data/
    tf_rules.yaml
    lecture_index.yaml
  tests/
    test_midterm_2025.py
    test_midterm_2024.py
    test_previous_midterms.py
  ui/
    cli.py
  notebooks/
    scratch.ipynb
  README.md
```

---

## 13. v1에서 반드시 구현해야 하는 기능

### 반드시 필요한 것
- T/F verdict + reason
- 3차 FRF 분해
- Lagrange EOM 생성
- equilibrium 찾기
- perturbation linearization
- symmetric EVP
- modal normalization
- proportional damping modal response
- markdown 풀이 출력

### 있으면 좋은 것
- 간단한 plot
- mode shape 시각화
- 입력 템플릿 자동 완성
- LaTeX 스타일 출력

### v1에서 굳이 안 해도 되는 것
- OCR
- 자유문장 전체 자동해석
- 범용 PDE solver
- FEM
- 웹 연동
- 대화형 LLM 의존 기능

---

## 14. 추천 개발 순서

### 1주차
- 프로젝트 뼈대 생성
- FRF 엔진 초안
- T/F rule DB 초안

### 2주차
- Lagrange engine 초안
- rotating hoop / triangle 템플릿 구현
- equilibrium / linearization 기능 구현

### 3주차
- modal engine 구현
- 2-DOF undamped 문제 회귀 테스트
- markdown report_engine 구현

### 4주차
- proportional damping 구현
- 2024 / 2025 midterm full regression
- CLI polishing

### 이후 확장
- state-space general solver
- gyroscopic solver
- nonconservative solver

---

## 15. 지금 시점의 최종 결론

현재 자료 기준으로 가장 효율적인 전략은 다음과 같다.

> **범용 수학 프로그램을 만들기보다, ME551 중간고사에 반복 출제되는 문제 유형을 직접 해결하는 템플릿형 수식 엔진을 만드는 것**

그리고 개발 우선순위는 아래가 맞다.

1. **T/F 개념 엔진**
2. **FRF 분해 엔진**
3. **Lagrange → 평형점 → 선형화 → 안정성 엔진**
4. **2-DOF modal analysis 엔진**
5. **proportional damping 엔진**
6. **gyroscopic / nonconservative / bi-orthogonality 확장**

즉 v1의 목표는

> **“시험 문제를 읽고 유형을 고른 뒤, 필요한 파라미터를 넣으면 답과 풀이가 바로 나오는 오프라인 solver”**

로 설정하는 것이 가장 현실적이다.

---

## 16. 참고 자료(이 plan을 세울 때 반영한 범위)

- Lecture 1: 과목 소개, ME551이 이론+계산 중심이며 중간고사가 Chap.1~4 중심이라는 구조
- Lecture 2~3: linear system, FRF, convolution, state equation, transition matrix
- Lecture 5~6: virtual work, D’Alembert, Hamilton, Lagrange
- Lecture 7~8: SDOF, structural damping, convolution, state equation of vibration
- Lecture 9~10: equilibrium, perturbation, linearization, stability
- Lecture 11: symmetric eigenvalue problem, orthogonality, modal analysis
- Lecture 12: gyroscopic systems
- Lecture 13: nonconservative systems, left/right eigenvectors, bi-orthogonality
- Midterm 2024 / 2025 / previous selected midterms: 반복되는 문제 유형 추출

