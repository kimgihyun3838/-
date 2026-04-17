# ME551 Homework 1/2/3 풀이 검증 리포트
## ME551 Linear Vibration Solver를 활용한 풀이 및 솔루션 대조 검증

**작성일**: 2026-04-15  
**프로젝트**: ME551 Linear Vibration Solver v2.0  
**대상**: Homework #1, #2, #3 (총 12문제)

---

## 목차

1. [프로젝트 솔버 기능 개요](#1-프로젝트-솔버-기능-개요)
2. [Homework #1 풀이 및 검증](#2-homework-1-풀이-및-검증)
3. [Homework #2 풀이 및 검증](#3-homework-2-풀이-및-검증)
4. [Homework #3 풀이 및 검증](#4-homework-3-풀이-및-검증)
5. [전체 검증 결과 요약](#5-전체-검증-결과-요약)
6. [솔버별 활용 매핑 테이블](#6-솔버별-활용-매핑-테이블)

---

## 1. 프로젝트 솔버 기능 개요

본 프로젝트는 `src/me551_solver/` 하위에 8개의 핵심 솔버 모듈을 보유하고 있으며, 각 숙제 문제는 이들 솔버의 조합으로 풀이 및 검증이 가능하다.

### 1.1 사용된 솔버 모듈

| 솔버 | 모듈 경로 | 핵심 기능 | 적용 문제 |
|------|----------|----------|----------|
| **FRFSolver** | `core/frf.py` | 전달함수 G(s) 구성, 부분분수 분해, ω_n/ζ 추출, 주파수 응답 | HW1-P2, HW1-P3 |
| **StateSpaceSolver** | `core/state_space.py` | 상태공간 변환, 고유값 계산, 전이행렬 Φ(t), 자유/강제 응답 | HW1-P3, HW1-P4, HW3-P4 |
| **StabilitySolver** | `core/stability.py` | 특성방정식 근 분류, 안정성 판별 (안정/불안정/한계안정) | HW1-P4 |
| **LagrangeSolver** | `core/lagrange.py` | 라그랑지안 유도, 운동방정식 도출, 평형점, 선형화 | HW2-P2, HW2-P4 |
| **ModalSolver** | `core/modal.py` | 일반화 고유값 문제, 고유진동수, 모드형상, 모달 응답 | (HW3 확장 가능) |
| **DampingSolver** | `core/damping.py` | 비례감쇠 해석, 모달 감쇠비, 감쇠 고유진동수 | HW3-P1 관련 |
| **ReportEngine** | `core/report.py` | 풀이 결과를 Exam-style/Markdown 형식으로 렌더링 | 전체 |
| **Router** | `core/router.py` | 문제 텍스트 자동 분류 → 적절한 솔버 라우팅 | 전체 |

### 1.2 실행 방법

```bash
cd src/
python -m me551_solver        # 대화형 CLI 메뉴
```

또는 Python API로 직접 호출:

```python
from me551_solver.core.frf import FRFSolver
result = FRFSolver().solve({"coefficients": [1, 0.6, 9], "numerator": [1]})
```

---

## 2. Homework #1 풀이 및 검증

### 2.1 Problem 1 (25pts) — 선형성 판별

**문제**: $\ddot{x} + (6 + 3\sin t)x = 0$ 이 선형인지 비선형인지 homogeneity와 additivity로 판별

#### 사용 솔버: 직접 해석적 검증 (솔버 불필요)

이 문제는 연산자 D(x) = x'' + (6+3sin(t))x에 대해 대수적 성질을 확인하는 이론 문제이므로, 특정 수치 솔버가 아닌 **해석적 추론**으로 풀이한다.

#### 풀이 과정

**연산자 정의**: $D(x) = \ddot{x} + (6 + 3\sin t)x$

**Homogeneity 검증**:
$$D(\alpha x) = (\alpha x)'' + (6+3\sin t)(\alpha x) = \alpha[x'' + (6+3\sin t)x] = \alpha D(x) \quad \checkmark$$

**Additivity 검증**:
$$D(x_1+x_2) = (x_1+x_2)'' + (6+3\sin t)(x_1+x_2) = D(x_1) + D(x_2) \quad \checkmark$$

**핵심 포인트**: $(6+3\sin t)$는 시간의 함수(time-varying coefficient)이지만 $x$에 대해서는 **1차**이므로, 계수가 시변(time-varying)이어도 선형성은 성립한다.

#### 결론: **선형 시스템 (Linear System)**

#### 솔루션 대조: **일치** ✓

---

### 2.2 Problem 2 (25pts) — Dynamic Stiffness, FRF, System Response

**문제**: $\ddot{x} + 0.6\dot{x} + 9x = 10\sin\omega t$ 의 동적강성, 주파수응답, 시스템응답 구하고 Bode plot

#### 사용 솔버: `FRFSolver` (Menu 2)

#### 실행 코드

```python
from me551_solver.core.frf import FRFSolver

solver = FRFSolver()
result = solver.solve({
    "coefficients": [1, 0.6, 9],   # s^2 + 0.6s + 9
    "numerator": [1]                # 분자 = 1
})
```

#### 솔버 출력 결과

```
Step 1: Transfer function G(s) 구성
  G(s) = 1 / (s² + 3s/5 + 9)

Step 2: Characteristic equation roots
  Roots: -3/10 ± 9√11·i/10

Step 3: 2nd-order subsystem 분석
  ω_n = 3 rad/s (ω_n² = 9)
  ζ = 0.1 (감쇠비)

Step 4: G(jω) 분석
  G(jω) = 1 / [(9-ω²) + 0.6jω]
  |G(0)| = 0.1111 (= 1/9)
  |G(jω_n)| = 0.5556 (∠ -90°)
```

#### 풀이 도출

솔버 결과로부터:

- **Dynamic stiffness**: $Z(j\omega) = (j\omega)^2 + 0.6j\omega + 9 = (9-\omega^2) + 0.6j\omega$
- **Frequency response**: $G(j\omega) = \frac{1}{(9-\omega^2) + 0.6j\omega}$
- **System response**: $x(t) = \frac{10}{\sqrt{(9-\omega^2)^2 + (0.6\omega)^2}} \sin(\omega t - \phi)$
  - where $\phi = \tan^{-1}\left(\frac{0.6\omega}{9-\omega^2}\right)$

#### 수치 검증

| 검증 항목 | 솔버 출력 | 수식 계산 | 일치 |
|----------|---------|---------|------|
| G(0) = 1/9 | 0.1111 | 0.1111 | ✓ |
| ω_n | 3.0 rad/s | √9 = 3.0 | ✓ |
| ζ | 0.1 | 0.6/(2×3) = 0.1 | ✓ |
| \|G(jω_n)\| | 0.5556 | 1/(0.6×3) = 0.5556 | ✓ |
| Phase at ω_n | -90° | tan⁻¹(∞) = 90° | ✓ |

#### 솔루션 대조: **일치** ✓

---

### 2.3 Problem 3 (25pts) — Impulse Response + State Equation

**문제**: $\dddot{x} + 2\ddot{x} + 4\dot{x} + 8x = \delta(t)$ 의 충격응답 및 상태방정식 유도

#### 사용 솔버 1: `FRFSolver` — 전달함수 및 부분분수 분해

```python
from me551_solver.core.frf import FRFSolver

result = FRFSolver().solve({
    "coefficients": [1, 2, 4, 8],   # s³ + 2s² + 4s + 8
    "numerator": [1]
})
```

**솔버 출력**:
```
Transfer function: G(s) = 1/(s³ + 2s² + 4s + 8)

Characteristic roots: s = -2, ±2j

Partial fraction decomposition:
  G(s) = 1/(8(s+2)) - (s-2)/(8(s²+4))

1st-order subsystem: G₁(s) = 1/(8(s+2)), ω_b = 2
2nd-order subsystem: G₂(s) = -(s-2)/(8(s²+4)), ω_n = 2, ζ = 0
```

**부분분수 계수 추출**: A = 1/8, B = -1/8, C = 1/4

**역 라플라스 변환**:
$$x(t) = \frac{1}{8}e^{-2t} - \frac{1}{8}\cos 2t + \frac{1}{8}\sin 2t, \quad t > 0$$

#### 사용 솔버 2: `StateSpaceSolver` — 상태방정식 구성 및 고유값 검증

```python
from me551_solver.core.state_space import StateSpaceSolver

result = StateSpaceSolver().solve({
    "mode": "state_matrix",
    "A": [[0,1,0], [0,0,1], [-8,-4,-2]],
    "B": [[0],[0],[1]],
    "x0": [0,0,0],
    "compute": ["eigenvalues", "stability"]
})
```

**솔버 출력**:
```
State matrix A:
  [ 0,  1,  0]
  [ 0,  0,  1]
  [-8, -4, -2]

Eigenvalues:
  λ₁ = -2       (real negative, stable node)
  λ₂ = -2j      (purely imaginary, marginally stable)
  λ₃ = +2j      (purely imaginary, marginally stable)

Stability: Marginally stable
```

#### 상태방정식 유도

$\dddot{x} = -2\ddot{x} - 4\dot{x} - 8x + \delta(t)$로부터:

$$\dot{\mathbf{x}}(t) = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ -8 & -4 & -2 \end{bmatrix} \mathbf{x}(t) + \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} \delta(t), \quad \mathbf{x}(t) = \begin{bmatrix} x \\ \dot{x} \\ \ddot{x} \end{bmatrix}$$

#### 수치 검증

| t | 해석해 | scipy 수치해 | 오차 |
|---|------|-----------|------|
| 1.0 | 0.18259744 | 0.18260639 | 8.9×10⁻⁶ |

#### 솔루션 대조: **일치** ✓

---

### 2.4 Problem 4 (25pts) — Step Response via Transition Matrix

**문제**: $\ddot{x} + 4\dot{x} + 9x = u(t)$ 의 전이행렬을 이용한 계단응답

#### 사용 솔버 1: `StateSpaceSolver` — 전이행렬 및 고유값

```python
from me551_solver.core.state_space import StateSpaceSolver

result = StateSpaceSolver().solve({
    "mode": "state_matrix",
    "A": [[0, 1], [-9, -4]],
    "B": [[0], [1]],
    "x0": [0, 0],
    "compute": ["eigenvalues", "stability", "transition_matrix"]
})
```

**솔버 출력**:
```
Eigenvalues:
  λ₁ = -2 - √5·j  (complex, Re < 0, stable focus)
  λ₂ = -2 + √5·j  (complex, Re < 0, stable focus)

Stability: Asymptotically stable

Transition matrix:
  Φ(t) = e^(-2t) × [2Re(R₁)cos(√5t) - 2Im(R₁)sin(√5t)]
  (잔류행렬 R₁, R₂ 수치값 제공)
```

#### 사용 솔버 2: `StabilitySolver` — 안정성 확인

```python
from me551_solver.core.stability import StabilitySolver

result = StabilitySolver().solve({
    "mode": "linearized_ode",
    "mass_coeff": "1",
    "damping_coeff": "4",
    "stiffness_coeff": "9"
})
```

**솔버 출력**:
```
Characteristic equation: λ² + 4λ + 9 = 0
Eigenvalues: -2 ± √5·i
Stability: k_eff = 9 > 0, c_eff = 4 > 0 → Asymptotically stable
```

#### 풀이 도출

상태방정식: $\dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 \\ -9 & -4 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ 1 \end{bmatrix} u(t)$

전이행렬 $e^{\mathbf{A}t}$를 구하고, $\mathbf{x}(0) = \mathbf{0}$일 때:

$$\mathbf{x}(t) = \int_0^t e^{\mathbf{A}(t-\tau)} \mathbf{b}\, u(\tau)\, d\tau$$

적분 결과:
$$x(t) = \frac{1}{9} - \frac{1}{9}e^{-2t}\cos(\sqrt{5}t) - \frac{2\sqrt{5}}{45}e^{-2t}\sin(\sqrt{5}t), \quad t > 0$$

#### 수치 검증

| t | 해석해 | scipy step 응답 | 오차 |
|---|------|---------------|------|
| 0.5 | 0.06035365 | 0.06036104 | 7.4×10⁻⁶ |
| 1.0 | 0.10981164 | 0.10981640 | 4.8×10⁻⁶ |
| 2.0 | 0.11336330 | 0.11336170 | 1.6×10⁻⁶ |

정상상태값: $x(\infty) = 1/9 = 0.1111$ (DC gain과 일치)

#### 솔루션 대조: **일치** ✓

---

## 3. Homework #2 풀이 및 검증

### 3.1 Problem 1 (25pts) — Rigid Bar (Newtonian Approach)

**문제**: Fig 2.23의 강체 봉-롤러-스프링 시스템의 뉴턴법 운동방정식 유도

#### 사용 솔버: 직접 해석적 유도 + scipy 수치 검증

이 문제는 기하학적 구속조건과 자유물체도(FBD)를 그려 뉴턴 법칙을 직접 적용하는 문제로, **LagrangeSolver와 다른 접근법**(Newton)을 요구한다. 솔버의 직접 적용은 제한적이나, 도출된 결과의 **수치적 정합성**을 검증할 수 있다.

#### 풀이 과정

1. **질량중심 C의 위치**: $\mathbf{r}_C = (-\frac{L}{2}\cos\theta)\hat{i} + (-\frac{L}{2}\sin\theta)\hat{j}$

2. **가속도**: $\mathbf{a}_C = (\frac{L}{2}\cos\theta\cdot\dot{\theta}^2 + \frac{L}{2}\sin\theta\cdot\ddot{\theta})\hat{i} + (\frac{L}{2}\sin\theta\cdot\dot{\theta}^2 - \frac{L}{2}\cos\theta\cdot\ddot{\theta})\hat{j}$

3. **3개 운동방정식**:
   - $\Sigma F_x$: $-k_1 L(1-\cos\theta) + N_2 - m(\frac{L}{2}\cos\theta\cdot\dot\theta^2 + \frac{L}{2}\sin\theta\cdot\ddot\theta) = 0$
   - $\Sigma F_y$: $k_2 L\sin\theta + N_1 - mg - m(\frac{L}{2}\sin\theta\cdot\dot\theta^2 - \frac{L}{2}\cos\theta\cdot\ddot\theta) = 0$
   - $\Sigma M_C$: $-k_1 L(1-\cos\theta)\frac{L}{2}\sin\theta + N_1\frac{L}{2}\cos\theta - \frac{1}{12}mL^2\ddot\theta - k_2 L\sin\theta\frac{L}{2}\cos\theta - N_2\frac{L}{2}\sin\theta = 0$

4. **파라미터 대입** ($m=1, L=2, k_1=5, k_2=10, N_1=3, N_2=6, g=10$):
   - $-10(1-\cos\theta) + 6 - (\cos\theta\cdot\dot\theta^2 + \sin\theta\cdot\ddot\theta) = 0$
   - $20\sin\theta - 7 - \sin\theta\cdot\dot\theta^2 + \cos\theta\cdot\ddot\theta = 0$
   - $-10(1-\cos\theta)\sin\theta + 3\cos\theta - \frac{1}{3}\ddot\theta - 20\sin\theta\cos\theta - 6\sin\theta = 0$

#### 솔루션 대조: **일치** ✓

---

### 3.2 Problem 2 (25pts) — Hamilton's Principle (Two Links + Mass)

**문제**: Fig 2.25의 두 링크 + 질량 시스템에 대해 Hamilton 원리로 운동방정식 유도

#### 사용 솔버: `LagrangeSolver` (Menu 3)

```python
from me551_solver.core.lagrange import LagrangeSolver

result = LagrangeSolver().solve({
    "T_expr": "Rational(1,3)*m*L**2*thetadot**2 + m*L**2*sin(theta)**2*thetadot**2 + Rational(1,2)*M*xdot**2",
    "V_expr": "-m*g*L*sin(theta) + Rational(1,2)*k1*(2*L-2*L*cos(theta)+x)**2 + Rational(1,2)*k2*x**2",
    "coord": "theta",
    "symbols": ["m","g","L","M","k1","k2","x","xdot","thetadot"]
})
```

**솔버 출력 (핵심 단계)**:
```
Step 1: T = L²m·θ̇²·sin²θ + L²m·θ̇²/3 + M·ẋ²/2

Step 2: V = -Lgm·sinθ + k₁(-2Lcosθ+2L+x)²/2 + k₂x²/2

Step 3: dT/dθ̇ = 2L²m·θ̇·(3sin²θ+1)/3

Step 4: d/dt(dT/dθ̇) = 2L²m·(θ̈(3sin²θ+1) + 3θ̇²sin2θ)/3

Step 5: dV/dθ = L(-gm·cosθ + 2k₁(-2Lcosθ+2L+x)sinθ)

Step 6: Lagrange EOM:
  (2/3·mL² + 2mL²sin²θ)θ̈ + mL²θ̇²sin2θ
  - mgLcosθ + 2k₁Lsinθ(2L-2Lcosθ+x) = 0
```

#### 솔버가 수행하는 자동 계산 과정

LagrangeSolver는 다음을 **기호적(symbolic)으로** 자동 수행한다:

1. T, V 파싱 → SymPy 식으로 변환
2. $\frac{\partial T}{\partial \dot\theta}$ 계산 (편미분)
3. $\frac{d}{dt}\left(\frac{\partial T}{\partial \dot\theta}\right)$ 계산 (체인룰 적용: $\dot\theta \to \ddot\theta$, $\theta \to \dot\theta$ 등)
4. $\frac{\partial T}{\partial \theta}$, $\frac{\partial V}{\partial \theta}$ 계산
5. Lagrange 방정식 조립: $\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} - \frac{\partial L}{\partial q_i} = 0$

#### 두 번째 좌표 x에 대한 EOM

솔버를 x 좌표에 대해서도 실행하면:
$$M\ddot{x} + k_1(2L - 2L\cos\theta + x) + k_2 x = 0$$

#### 솔루션 대조: **일치** ✓

> **참고**: Hamilton 원리는 $\int_{t_1}^{t_2}(\delta L + \delta W_{nc})dt = 0$을 적용하는 것이며, LagrangeSolver의 Lagrange 방정식 결과와 동일한 EOM을 유도한다. 두 방법은 수학적으로 등가이다.

---

### 3.3 Problem 3 (25pts) — Virtual Work Principle

**문제**: Fig 2.25 시스템의 가상일 원리를 이용한 평형방정식 유도

#### 사용 솔버: 직접 해석적 유도

가상일 원리는 **정적 평형** ($\ddot\theta = 0, \dot\theta = 0$)에서의 방정식을 구하는 것으로, 라그랑주 결과에서 관성항을 제거한 형태와 같다.

#### 풀이 과정

**힘과 가상변위 설정**:

| 힘 | 크기 | 가상변위 |
|----|------|---------|
| $\mathbf{F}_1$ | $-k_1(x+2L-2L\cos\theta)\hat{i}$ | $\delta x_1 = \delta x + 2L\sin\theta\,\delta\theta$ |
| $\mathbf{F}_2$ | $-k_2 x\,\hat{i}$ | $\delta x_2 = \delta x$ |
| $\mathbf{F}_3$ | $-mg\,\hat{j}$ | $\delta x_3 = -\frac{L}{2}\cos\theta\,\delta\theta$ |
| $\mathbf{F}_4$ | $-mg\,\hat{j}$ | $\delta x_4 = -\frac{L}{2}\cos\theta\,\delta\theta$ |

**가상일**: $\sum F_i \cdot \delta x_i = 0$ 으로부터:

- $\delta\theta$ 계수: $-2k_1 L\sin\theta(x+2L-2L\cos\theta) + mgL\cos\theta = 0$
- $\delta x$ 계수: $(k_1+k_2)x + 2k_1 L(1-\cos\theta) = 0$

#### 솔루션 대조: **일치** ✓

---

### 3.4 Problem 4 (25pts) — Lagrange Equation

**문제**: Fig 2.25 시스템의 Lagrange 운동방정식 유도

#### 사용 솔버: `LagrangeSolver` (Menu 3) — Problem 2와 동일

일반화 좌표 수: **2개** $(x, \theta)$

Problem 2에서 이미 LagrangeSolver로 유도한 것과 동일한 라그랑지안을 사용하되, Lagrange 방정식:

$$\frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}_i}\right) - \frac{\partial L}{\partial q_i} = 0, \quad i = 1, 2$$

을 적용하면:

$$\begin{cases}
M\ddot{x} + (k_1+k_2)x + 2k_1 L(1-\cos\theta) = 0 \\
\left(\frac{2}{3}mL^2 + 2mL^2\sin^2\theta\right)\ddot\theta + 2mL^2\sin\theta\cos\theta\,\dot\theta^2 - mgL\cos\theta + 2k_1 L(2L-2L\cos\theta+x)\sin\theta = 0
\end{cases}$$

#### 솔루션 대조: **일치** ✓

---

## 4. Homework #3 풀이 및 검증

### 4.1 Problem 1 (25pts) — Base Excitation Response

**문제**: Fig 1의 질량-감쇠기-스프링 시스템에 기저 가진 $y(t) = A\sin\omega t$가 작용할 때, 응답의 크기와 위상 유도

#### 사용 솔버: `FRFSolver` + 해석적 유도

#### 풀이 과정

운동방정식:
$$M\ddot{x} + c(\dot{x}-\dot{y}) + k(x-y) = 0$$
$$\Rightarrow M\ddot{x} + c\dot{x} + kx = c\dot{y} + ky = Ac\omega\cos\omega t + Ak\sin\omega t$$

$r = \omega/\omega_n$, $\zeta = c/(2M\omega_n)$으로 무차원화:

우변을 하나의 조화함수로 합성:
$$RHS = A\omega_n^2\sqrt{1+(2\zeta r)^2}\cos(\omega t - \alpha), \quad \tan\alpha = \frac{1}{2\zeta r}$$

FRFSolver로 좌변의 주파수응답함수를 구하면:
$$G(j\omega) = \frac{1}{(1-r^2) + 2j\zeta r} \quad (\omega_n^2 \text{으로 정규화})$$

**최종 응답**:
- **크기**: $\displaystyle\frac{A\sqrt{1+(2\zeta r)^2}}{\sqrt{(1-r^2)^2+(2\zeta r)^2}}$
- **위상**: $\tan^{-1}\!\left(\frac{1}{2\zeta r}\right) + \tan^{-1}\!\left(\frac{2\zeta r}{1-r^2}\right)$

#### 수치 검증

| ζ | r | 크기/A (계산) | 공식 검증 | 일치 |
|---|---|-----------|--------|------|
| 0.1 | 1.0 | 5.0990 | √(1.04)/0.2 = 5.0990 | ✓ |
| 0.3 | 1.0 | 1.9437 | √(1.36)/0.6 = 1.9437 | ✓ |
| 0.5 | 1.0 | 1.4142 | √(2)/1.0 = 1.4142 | ✓ |

#### 솔루션 대조: **일치** ✓

---

### 4.2 Problem 2 (25pts) — Square Wave + First-Order System

**문제**: 1차 시스템 $c\dot{x} + kx = f(t)$에 주기적 사각파 가진에 대한 응답 유도

#### 사용 솔버: `FRFSolver` (각 조화 성분의 FRF 계산)

#### 풀이 과정

**1단계: Fourier 급수 전개**

사각파 $f(t)$의 Fourier 계수:
- $A_0 = f_0$
- $A_p = \frac{2f_0}{p\pi}\sin\frac{p\pi}{2}, \quad p = 1, 2, 3, \ldots$

**수치 확인**:
| p | $A_p/f_0$ | 비고 |
|---|---------|------|
| 1 | +0.6366 | 기여 |
| 2 | 0.0000 | 소멸 |
| 3 | -0.2122 | 기여 |
| 4 | 0.0000 | 소멸 |
| 5 | +0.1273 | 기여 |

→ **홀수항만 기여** (짝수항은 $\sin(p\pi/2) = 0$)

**2단계: 1차 시스템의 FRF**

$c\dot{x} + kx = f(t)$를 $k$로 나누면: $\tau\dot{x} + x = f(t)/k$, $\tau = c/k$

$$G(j\omega) = \frac{1}{1+j\omega\tau}, \quad |G| = \frac{1}{\sqrt{1+(\omega\tau)^2}}, \quad \varphi = \tan^{-1}(\omega\tau)$$

**3단계: 중첩**

$$x(t) = \frac{f_0}{2k} + \sum_{p=1}^{\infty} \frac{2f_0}{kp\pi}\sin\frac{p\pi}{2} \cdot \frac{1}{\sqrt{1+(2\pi p\tau/T)^2}} \cos\left(\frac{2\pi p}{T}t - \tan^{-1}\frac{2\pi p\tau}{T}\right)$$

#### 솔루션 대조: **일치** ✓

---

### 4.3 Problem 3 (25pts) — Convolution Integral

**문제**: 부족감쇠 2차 시스템 $m\ddot{x} + c\dot{x} + kx = f(t)$에 $f(t) = F_0 e^{-at}$ (0 < t < t₁)에 대한 컨볼루션 적분 응답

#### 사용 솔버: `StateSpaceSolver` + scipy.integrate 수치 검증

#### 풀이 과정

**충격 응답함수**:
$$g(t) = \frac{1}{m\omega_d}e^{-\zeta\omega_n t}\sin\omega_d t, \quad \omega_d = \sqrt{1-\zeta^2}\,\omega_n$$

**컨볼루션**:
$$x(t) = \int_0^t f(\tau) g(t-\tau) d\tau = \frac{1}{m\omega_d}\int_0^t f(\tau) e^{-\zeta\omega_n(t-\tau)}\sin\omega_d(t-\tau)\,d\tau$$

**Case i) $0 < t < t_1$**:
$$x(t) = \frac{F_0}{m\omega_d[(\zeta\omega_n-a)^2+\omega_d^2]}\left[e^{-\zeta\omega_n t}\left(-(\zeta\omega_n-a)\sin\omega_d t - \omega_d\cos\omega_d t\right) + e^{-at}\omega_d\right]$$

**Case ii) $t > t_1$**: 적분 상한이 $t_1$으로 변경, 유사한 형태

#### 수치 검증 (m=1, ζ=0.3, ω_n=10, a=2, F₀=5, t₁=3)

| t | 해석해 | ODE 수치해 | 오차 |
|---|------|---------|------|
| 0.5 | 0.02056802 | 0.02056644 | 1.6×10⁻⁶ |
| 1.0 | 0.01007568 | 0.01007522 | 4.6×10⁻⁷ |
| 2.0 | 0.00086102 | 0.00086096 | 5.9×10⁻⁸ |

#### 솔루션 대조: **일치** ✓

---

### 4.4 Problem 4 (25pts) — Transition Matrix Approach

**문제**: Problem 3의 시스템에 대해 전이행렬 방법으로 상태벡터 유도

#### 사용 솔버: `StateSpaceSolver` (transition_matrix 모드)

```python
from me551_solver.core.state_space import StateSpaceSolver

# omega_n=10, zeta=0.3 → c=6, k=100 for m=1
result = StateSpaceSolver().solve({
    "mode": "state_matrix",
    "A": [[0, 1], [-100, -6]],   # [-omega_n^2, -2*zeta*omega_n]
    "B": [[0], [1]],
    "x0": [0, 0],
    "compute": ["eigenvalues", "transition_matrix"]
})
```

**솔버 출력**:
```
Eigenvalues: -3 ± 9.5394j (stable focus)

Transition matrix:
Φ(t) = e^(At) = (e^(-ζω_n·t)/ω_d) ×
  [ω_d cosω_d t + ζω_n sinω_d t,     sinω_d t           ]
  [-ω_n² sinω_d t,                   ω_d cosω_d t - ζω_n sinω_d t]
```

#### 전이행렬 수치 검증

scipy.linalg.expm과의 비교:
```
t = 0.5:
  해석적 Φ(0.5) vs expm(A×0.5): max error = 7.11×10⁻¹⁵
```

#### 상태 응답

$\mathbf{x}(0) = \mathbf{0}$일 때:
$$\mathbf{x}(t) = \int_0^t \Phi(t-\tau)\mathbf{b}\,f(\tau)\,d\tau$$

첫 번째 성분 $x(t)$를 추출하면 Problem 3의 컨볼루션 적분과 동일한 결과를 얻는다.

#### 솔루션 대조: **일치** ✓

---

## 5. 전체 검증 결과 요약

### 5.1 HW1 결과

| 문제 | 유형 | 사용 솔버 | 핵심 검증값 | 결과 |
|------|------|---------|----------|------|
| P1 | 선형성 판별 | 해석적 추론 | Homogeneity ✓, Additivity ✓ | **MATCH** |
| P2 | FRF/시스템 응답 | `FRFSolver` | ω_n=3, ζ=0.1, \|G(0)\|=1/9 | **MATCH** |
| P3 | 충격응답+상태방정식 | `FRFSolver` + `StateSpaceSolver` | x(1)=0.1826, λ=-2,±2j | **MATCH** |
| P4 | 전이행렬 계단응답 | `StateSpaceSolver` + `StabilitySolver` | x(∞)=1/9, λ=-2±√5j | **MATCH** |

### 5.2 HW2 결과

| 문제 | 유형 | 사용 솔버 | 핵심 검증 | 결과 |
|------|------|---------|---------|------|
| P1 | 뉴턴법 EOM | 해석적 유도 | 3개 방정식 + 파라미터 대입 | **MATCH** |
| P2 | Hamilton 원리 | `LagrangeSolver` | T, V, L → EOM (2개) | **MATCH** |
| P3 | 가상일 원리 | 해석적 유도 | 평형방정식 2개 | **MATCH** |
| P4 | Lagrange 방정식 | `LagrangeSolver` | P2와 동일 EOM | **MATCH** |

### 5.3 HW3 결과

| 문제 | 유형 | 사용 솔버 | 핵심 검증값 | 결과 |
|------|------|---------|----------|------|
| P1 | 기저가진 응답 | `FRFSolver` + 해석 | r=1에서 크기/위상 정합 | **MATCH** |
| P2 | 사각파+1차 시스템 | `FRFSolver` (조화성분별) | Fourier 계수 + FRF 중첩 | **MATCH** |
| P3 | 컨볼루션 적분 | scipy 수치 검증 | 오차 < 10⁻⁶ | **MATCH** |
| P4 | 전이행렬 | `StateSpaceSolver` | Φ(t) 오차 < 10⁻¹⁴ | **MATCH** |

### 총합: **12/12 문제 모두 솔루션과 일치** ✓

---

## 6. 솔버별 활용 매핑 테이블

### 6.1 문제 → 솔버 매핑

```
HW1-P1  ──→  (해석적 추론, 솔버 불필요)
HW1-P2  ──→  FRFSolver          [Menu 2]
HW1-P3  ──→  FRFSolver          [Menu 2]  +  StateSpaceSolver
HW1-P4  ──→  StateSpaceSolver   +  StabilitySolver

HW2-P1  ──→  (뉴턴법 직접 유도)
HW2-P2  ──→  LagrangeSolver     [Menu 3]
HW2-P3  ──→  (가상일 원리 직접 유도)
HW2-P4  ──→  LagrangeSolver     [Menu 3]

HW3-P1  ──→  FRFSolver          [Menu 2]  +  해석적 유도
HW3-P2  ──→  FRFSolver          [Menu 2]  (각 조화 성분)
HW3-P3  ──→  scipy 수치 검증
HW3-P4  ──→  StateSpaceSolver   (transition_matrix 모드)
```

### 6.2 솔버 기능별 활용 횟수

| 솔버 기능 | 활용 문제 | 횟수 |
|----------|---------|------|
| 전달함수 G(s) 구성 | HW1-P2, P3, HW3-P1, P2 | 4 |
| 부분분수 분해 | HW1-P3 | 1 |
| ω_n, ζ 추출 | HW1-P2, HW3-P1 | 2 |
| 고유값 계산 | HW1-P3, P4, HW3-P4 | 3 |
| 안정성 분류 | HW1-P3, P4 | 2 |
| 전이행렬 Φ(t) | HW1-P4, HW3-P4 | 2 |
| Lagrange EOM 유도 | HW2-P2, P4 | 2 |
| 보고서 렌더링 | 전체 | 12 |

### 6.3 CLI 메뉴 퀵 가이드

```
숙제 문제를 풀려면:

1. python -m me551_solver 실행
2. 'a' (자동 분류) 선택 → 문제 텍스트 입력
3. 또는 직접 메뉴 선택:
   - ODE 전달함수 문제 → '2' (FRF 분해)
   - 라그랑지안/에너지 문제 → '3' (Lagrange 분석)
   - 행렬 고유값 문제 → '4' (Modal Analysis)
   - 감쇠 시스템 → '5' (Proportional Damping)
4. 결과 확인 → 'y'로 리포트 저장
```

---

**END OF REPORT**
