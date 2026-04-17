# ME551 2025 Midterm Exam Solutions

> Solved using ME551 Solver Project
> Date: 2026-04-13

---

## Problem 1 (30 pts) — True/False

### (a) Linear System Analysis (5pts)

> Response to arbitrary excitation can be obtained by either convolution integral of
> impulse response function and input in the time domain; or multiplication of FRF with
> input in the frequency domain in both cases of free vibration and forced vibration.

**Answer: FALSE**

Convolution integral과 FRF 곱셈은 **강제진동(forced vibration)에만** 적용된다.
자유진동(free vibration)은 외부 입력 f(t)가 없으므로 convolution/FRF가 아닌
초기조건(initial conditions)으로부터 응답을 구해야 한다.

- 강제진동: x(t) = h(t) * f(t), X(w) = H(w)F(w)
- 자유진동: 초기조건 x(0), x'(0)에 의한 응답

---

### (b) State Equation (5pts)

> A second-order linear MIMO system with N-dof can be formulated into
> 1st order N-dimensional state equation.

**Answer: FALSE**

2차 N-DOF 시스템은 **2N차원** 1차 상태방정식으로 변환된다.
각 DOF마다 변위와 속도 2개의 상태변수가 필요하므로:
- 상태벡터: x = [q; q_dot]^T (2N x 1)
- 상태행렬: A (2N x 2N)

**솔버 검증** (tf_rules.yaml rule_011):
> "A p-th order MIMO system with n inputs/outputs converted to first-order state form
> has a state vector of dimension p*n." → 2nd order, N-dof → dimension = 2N

---

### (c) Principle of Virtual Work (5pts)

> In the Principle of Virtual Work, effects of internal and constraint forces are eliminated.

**Answer: TRUE**

이상적인 구속(frictionless, rigid links 등)의 구속력은 가상변위에 수직이므로
가상일이 0이 되어 자동으로 소거된다. 내력(action-reaction pair)도 뉴턴 제3법칙에
의해 서로 상쇄되어 가상일 합이 0이다.

---

### (d) Structural Damping (5pts)

> Structural damping yields the same energy loss per cycle with equivalent damping.
> It can be used in either transient or harmonic analysis.

**Answer: FALSE**

Structural damping은 **조화 가진과 조화 응답에서만** 유효하다.
과도(transient) 해석에서는 사용할 수 없다.

등가 점성감쇠로 변환할 때 에너지 소산이 주파수에 독립적인 성질은
조화 운동을 전제로 유도된 것이다.

---

### (e) Eigenvalue Problem (5pts)

> For a real symmetric positive definite matrix A, we have real positive eigenvalues,
> and eigenvectors are also real-valued and orthogonal to each other.

**Answer: TRUE**

실수 대칭 양정치(SPD) 행렬의 성질:
- 고유값: 모두 실수이고 양수
- 고유벡터: 모두 실수이고 서로 직교 (spectral theorem)
- A = Q Λ Q^T (Q: 직교행렬, Λ: 양의 대각행렬)

---

### (f) Modal Analysis (5pts)

> Gyroscopic asymmetric system cannot be diagonalized using orthogonality of its
> right and left eigenvectors in modal analysis.

**Answer: FALSE**

자이로스코픽 비대칭 시스템은 좌고유벡터(left eigenvectors)와
우고유벡터(right eigenvectors)의 **쌍직교성(bi-orthogonality)**을 이용하여
대각화할 수 있다.

대칭 시스템: U^T M U = I (한 종류의 고유벡터)
비대칭 시스템: Ψ^T Φ (좌/우 고유벡터 쌍으로 대각화)

---

## Problem 2 (10 pts) — FRF Decomposition

> x''' + 3x'' + 6x' + 8x = f(t)

### (a) Decompose G₃(jw) into G₁(jw) + G₂(jw) (5pts)

**솔버 사용**: `FRFDecomposeSolver({'ode_string': "x''' + 3x'' + 6x' + 8x = f(t)"})`

#### Step 1: Transfer Function

$$G(s) = \frac{1}{s^3 + 3s^2 + 6s + 8}$$

#### Step 2: Characteristic Roots

$$s^3 + 3s^2 + 6s + 8 = 0$$
$$\text{Roots: } s = -2, \quad s = -\frac{1}{2} \pm \frac{\sqrt{15}}{2}j$$

#### Step 3: Factorization

$$s^3 + 3s^2 + 6s + 8 = (s+2)(s^2+s+4)$$

#### Step 4: Partial Fraction Decomposition

$$G(s) = \frac{1}{(s+2)(s^2+s+4)} = \underbrace{\frac{1}{6(s+2)}}_{G_1(s)} + \underbrace{\frac{-(s-1)}{6(s^2+s+4)}}_{G_2(s)}$$

#### Step 5: jw Domain

**G₁(jw):**

$$G_1(j\omega) = \frac{1}{(12) + j(6\omega)}$$

$$\text{Re}[G_1] = \frac{1}{3(\omega^2 + 4)}, \quad \text{Im}[G_1] = \frac{-\omega}{6\omega^2 + 24}$$

**G₂(jw):**

$$G_2(j\omega) = \frac{(1) + j(-\omega)}{(24 - 6\omega^2) + j(6\omega)}$$

$$\text{Re}[G_2] = \frac{2 - \omega^2}{3(\omega^4 - 7\omega^2 + 16)}, \quad \text{Im}[G_2] = \frac{\omega(\omega^2 - 5)}{6(\omega^4 - 7\omega^2 + 16)}$$

#### Characteristic Values

| Function | |G(0)| | |G(jw_b)| | |G(jw_n)| | |G(j∞)| |
|----------|--------|-----------|-----------|---------|
| G (total) | 0.125 | 0.1768 (∠-135°) | 0.1768 (∠-135°) | 0 |
| G₁ | 0.0833 | 0.0589 (∠-45°) | 0.0589 (∠-45°) | 0 |
| G₂ | 0.0417 | 0.1863 (∠-153°) | 0.1863 (∠-153°) | 0 |

*Note: 이 시스템에서 w_b = w_n = 2 rad/s이므로 두 열의 값이 동일*

### (b) Find w_b, w_n, zeta (5pts)

| Parameter | Value |
|-----------|-------|
| **w_b (break frequency)** | **2 rad/s** |
| **w_n (natural frequency)** | **2 rad/s** |
| **zeta (damping ratio)** | **0.25** |

- 1차 서브시스템: G₁(s) = 1/(6(s+2)) → w_b = 2 rad/s, tau = 0.5 s
- 2차 서브시스템: G₂(s) = -(s-1)/(6(s²+s+4)) → w_n = sqrt(4) = 2 rad/s, zeta = 1/(2*sqrt(4)) = 1/4 = 0.25

---

## Problem 3 (25 pts) — Lagrange Equation

> Rotating triangle (height L, width 2L), angular speed Omega, bead mass m, coordinate r.

### (a) Derive EOM using Lagrange equation (10pts)

**솔버 사용**: `LagrangeSolver` with T = m/2*(rdot² + Omega²r²), V = mgr

#### Kinetic Energy

$$T = \frac{1}{2}m(\dot{r}^2 + \Omega^2 r^2)$$

- 반경방향 속도: r_dot
- 접선방향 속도: Omega * r (회전에 의한)

#### Potential Energy

$$V = mgr$$

- 비드가 삼각형 위로 갈수록 높이가 증가 (기하학적 관계에서 V = mgr)

#### Lagrange Equation Derivation

$$\frac{\partial T}{\partial \dot{r}} = m\dot{r}$$

$$\frac{d}{dt}\frac{\partial T}{\partial \dot{r}} = m\ddot{r}$$

$$\frac{\partial T}{\partial r} = m\Omega^2 r$$

$$\frac{\partial V}{\partial r} = mg$$

$$\frac{d}{dt}\frac{\partial T}{\partial \dot{r}} - \frac{\partial T}{\partial r} + \frac{\partial V}{\partial r} = 0$$

$$\boxed{m\ddot{r} - m\Omega^2 r + mg = 0}$$

또는: $\ddot{r} - \Omega^2 r + g = 0$

### (b) Equilibrium position (5pts)

정적 평형: $\dot{r} = 0, \ddot{r} = 0$

$$-\Omega^2 r + g = 0$$

$$\boxed{r_{eq} = \frac{g}{\Omega^2}}$$

평형점이 삼각형 내부에 존재하려면 $r_{eq} < L$:

$$\frac{g}{\Omega^2} < L \quad \Rightarrow \quad \boxed{\Omega^2 > \frac{g}{L}}$$

### (c) Linearized EOM (5pts)

$r = r_{eq} + \delta r$ 대입, 선형화:

$$m\delta\ddot{r} + \underbrace{(-m\Omega^2)}_{k_{eff}} \delta r = 0$$

$$\boxed{\delta\ddot{r} - \Omega^2 \delta r = 0}$$

- m_eff = m
- c_eff = 0 (감쇠 없음)
- k_eff = -mOmega² (음수!)

### (d) Stability (5pts)

$$\boxed{\text{UNSTABLE (divergence)}}$$

**이유**: 유효 강성 k_eff = -mOmega² < 0 (항상 음수)

특성방정식: lambda² - Omega² = 0 → lambda = ±Omega

- lambda = +Omega > 0 인 양의 실수 고유값이 존재
- 이는 divergence 불안정 (비진동적 발산)
- 물리적 해석: 원심력이 복원력보다 항상 크므로 평형점에서 벗어나면 발산

---

## Problem 4 (35 pts) — Modal Analysis + Proportional Damping

> m₁ = m₂ = 9 kg, k₁ = k₂ = k₃ = 36 N/m
> IC: q₁(0) = q₂(0) = 1 m, q̇₁(0) = q̇₂(0) = 0 m/s

### (a) Matrix-vector EOM & eigenvalue problem (5pts)

**솔버 사용**: `ModalSolver` with M = [[9,0],[0,9]], K = [[72,-36],[-36,72]]

#### Mass & Stiffness Matrices

$$M = \begin{bmatrix} 9 & 0 \\ 0 & 9 \end{bmatrix}, \quad K = \begin{bmatrix} k_1+k_2 & -k_2 \\ -k_2 & k_2+k_3 \end{bmatrix} = \begin{bmatrix} 72 & -36 \\ -36 & 72 \end{bmatrix}$$

$$M\ddot{q} + Kq = 0$$

Harmonic solution $q(t) = u e^{j\omega t}$ 대입:

$$\boxed{Ku = \lambda Mu, \quad \lambda = \omega^2}$$

### (b) Standard eigenvalue problem (5pts)

Coordinate transform: $u = M^{-1/2} v$

$$A = M^{-1/2} K M^{-1/2} = \begin{bmatrix} 8 & -4 \\ -4 & 8 \end{bmatrix}$$

$$Av = \lambda v$$

**Eigenvalues:**

$$\boxed{\lambda_1 = 4, \quad \lambda_2 = 12}$$

$$\omega_1 = \sqrt{4} = 2 \text{ rad/s}, \quad \omega_2 = \sqrt{12} = 2\sqrt{3} \approx 3.464 \text{ rad/s}$$

**Eigenvectors (of A):**

$$v_1 = \begin{bmatrix} 0.7071 \\ 0.7071 \end{bmatrix}, \quad v_2 = \begin{bmatrix} 0.7071 \\ -0.7071 \end{bmatrix}$$

### (c) Mass-normalized modeshapes & orthogonality (5pts)

$$u_i = M^{-1/2} v_i$$

$$U = [u_1, u_2] = \begin{bmatrix} 0.2357 & 0.2357 \\ 0.2357 & -0.2357 \end{bmatrix}$$

**Orthogonality check:**

$$U^T M U = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} = I \quad \checkmark$$

$$U^T K U = \begin{bmatrix} 4 & 0 \\ 0 & 12 \end{bmatrix} = \text{diag}(\lambda_1, \lambda_2) \quad \checkmark$$

### (d) Proportional damping — decoupled modal equations (5pts)

**솔버 사용**: `DampingSolver` with alpha=0.2, beta=0.1

$$C = 0.2M + 0.1K = \begin{bmatrix} 9.0 & -3.6 \\ -3.6 & 9.0 \end{bmatrix}$$

Modal expansion $q(t) = u_1 \eta_1(t) + u_2 \eta_2(t)$, 양변에 $U^T$ 곱:

$$U^T M U \ddot{\eta} + U^T C U \dot{\eta} + U^T K U \eta = 0$$

$$I \ddot{\eta} + \text{diag}(2\zeta_1\omega_1, 2\zeta_2\omega_2) \dot{\eta} + \text{diag}(\omega_1^2, \omega_2^2) \eta = 0$$

Modal damping ratios: $2\zeta_r \omega_r = \alpha + \beta \omega_r^2$

| Mode | $\omega_r$ | $\alpha + \beta\omega_r^2$ | $\zeta_r$ | $\omega_{d,r}$ |
|------|-----------|------------------------|----------|--------------|
| 1 | 2 | 0.2 + 0.1(4) = 0.6 | **0.15** | 1.9774 |
| 2 | 3.464 | 0.2 + 0.1(12) = 1.4 | **0.2021** | 3.3926 |

**Decoupled equations:**

$$\boxed{\ddot{\eta}_1 + 0.6\dot{\eta}_1 + 4\eta_1 = 0}$$
$$\boxed{\ddot{\eta}_2 + 1.4\dot{\eta}_2 + 12\eta_2 = 0}$$

### (e) Modal response (10pts)

**Modal initial conditions:**

$$\eta(0) = U^T M q(0) = \begin{bmatrix} 0.2357 & 0.2357 \\ 0.2357 & -0.2357 \end{bmatrix} \begin{bmatrix} 9 & 0 \\ 0 & 9 \end{bmatrix} \begin{bmatrix} 1 \\ 1 \end{bmatrix} = \begin{bmatrix} 4.2426 \\ 0 \end{bmatrix}$$

$$\dot{\eta}(0) = U^T M \dot{q}(0) = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$

**Mode 1** ($\eta_1(0) = 4.2426$, $\dot{\eta}_1(0) = 0$, underdamped):

$$\eta_1(t) = e^{-\sigma_1 t}\left[A_1 \cos(\omega_{d,1} t) + B_1 \sin(\omega_{d,1} t)\right]$$

where $\sigma_1 = \zeta_1 \omega_1 = 0.3$

$$A_1 = \eta_1(0) = 4.2426$$
$$B_1 = \frac{\dot{\eta}_1(0) + \sigma_1 \eta_1(0)}{\omega_{d,1}} = \frac{0 + 0.3 \times 4.2426}{1.9774} = 0.6437$$

$$\boxed{\eta_1(t) = e^{-0.3t}\left[4.2426\cos(1.9774t) + 0.6437\sin(1.9774t)\right]}$$

**Mode 2** ($\eta_2(0) = 0$, $\dot{\eta}_2(0) = 0$):

$$\boxed{\eta_2(t) = 0}$$

(초기조건이 모두 0이므로 2번 모드는 여기되지 않음)

### (f) Physical response (5pts)

$$q(t) = U\eta(t) = u_1 \eta_1(t) + u_2 \eta_2(t)$$

$\eta_2(t) = 0$이므로:

$$q_1(t) = u_{11}\eta_1(t) = 0.2357 \times \eta_1(t)$$
$$q_2(t) = u_{21}\eta_1(t) = 0.2357 \times \eta_1(t)$$

$$\boxed{q_1(t) = q_2(t) = e^{-0.3t}\left[\cos(1.9774t) + 0.1518\sin(1.9774t)\right]}$$

**물리적 해석**: 초기조건 q₁(0) = q₂(0) = 1 (동위상)이므로 1번 모드(동위상 모드)만 여기되고,
2번 모드(역위상 모드)는 여기되지 않는다. 따라서 q₁(t) = q₂(t)로 두 질량이 항상 같은 운동을 한다.

---

## Summary of Answers

### Problem 1
| | Answer |
|---|--------|
| (a) | **FALSE** — convolution/FRF는 강제진동에만 적용 |
| (b) | **FALSE** — 2N차원 (N이 아님) |
| (c) | **TRUE** — 구속력과 내력은 가상일에서 소거됨 |
| (d) | **FALSE** — 조화 해석에서만 유효 |
| (e) | **TRUE** — SPD 행렬의 성질 |
| (f) | **FALSE** — 좌/우 고유벡터 쌍직교성으로 대각화 가능 |

### Problem 2
| Parameter | Value |
|-----------|-------|
| w_b | 2 rad/s |
| w_n | 2 rad/s |
| zeta | 0.25 |
| G₁(s) | 1/(6(s+2)) |
| G₂(s) | -(s-1)/(6(s²+s+4)) |

### Problem 3
| Item | Result |
|------|--------|
| EOM | m*r'' - m*Omega²*r + m*g = 0 |
| r_eq | g/Omega² |
| 조건 | Omega² > g/L |
| 선형화 | delta_r'' - Omega²*delta_r = 0 |
| 안정성 | **UNSTABLE** (k_eff = -mOmega² < 0, divergence) |

### Problem 4
| Item | Value |
|------|-------|
| lambda₁, lambda₂ | 4, 12 |
| omega₁, omega₂ | 2, 2sqrt(3) rad/s |
| U | [[0.2357, 0.2357], [0.2357, -0.2357]] |
| zeta₁, zeta₂ | 0.15, 0.2021 |
| omega_d1, omega_d2 | 1.9774, 3.3926 rad/s |
| eta₁(t) | e^(-0.3t)[4.2426cos(1.9774t) + 0.6437sin(1.9774t)] |
| eta₂(t) | 0 |
| q₁(t) = q₂(t) | e^(-0.3t)[cos(1.9774t) + 0.1518sin(1.9774t)] |
