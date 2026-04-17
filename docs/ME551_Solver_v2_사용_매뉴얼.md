# ME551 Solver v2.0 — 사용 매뉴얼

> 최종 업데이트: 2026-04-13
> 과목: ME551 Linear Vibration Engineering (선형진동공학)

---

## 목차

1. [시작하기](#1-시작하기)
2. [메뉴 a — 자동 분류](#2-메뉴-a--자동-분류-문제-유형-자동-판별)
3. [메뉴 1 — T/F 개념 판별](#3-메뉴-1--tf-개념-판별-v2)
4. [메뉴 2 — FRF 분해 (전달함수)](#4-메뉴-2--frf-분해-전달함수)
5. [메뉴 3 — Lagrange 분석 (1-DOF / N-DOF)](#5-메뉴-3--lagrange-분석)
6. [메뉴 4 — N-DOF Modal Analysis](#6-메뉴-4--n-dof-modal-analysis)
7. [메뉴 5 — Proportional Damping Response](#7-메뉴-5--proportional-damping-response)
8. [메뉴 6 — 확장 (Gyroscopic / Nonconservative)](#8-메뉴-6--확장-gyroscopic--nonconservative)
9. [메뉴 7~9 — 빈출 템플릿](#9-메뉴-79--빈출-템플릿)
10. [결과 저장 및 그래프](#10-결과-저장-및-그래프)
11. [v2 변경사항 요약](#11-v2-변경사항-요약)
12. [트러블슈팅](#12-트러블슈팅)

---

## 1. 시작하기

### 1.1 실행 방법

```bash
cd C:/Users/김기현/Desktop/선형진동공학/src
python -m me551_solver
```

### 1.2 메인 메뉴

```
=============================================
   ME551 Midterm Solver v2.0
   Linear Vibration Engineering
=============================================

------ Main Menu ----------------------------
 a. 자동 분류 (문제 입력 -> 유형 판별 -> 풀이)
 ---
 1. T/F 개념 판별
 2. FRF 분해 (전달함수)
 3. Lagrange -> 평형점 -> 선형화 -> 안정성
 4. N-DOF Modal Analysis
 5. Proportional Damping Response
 6. 확장 (Gyroscopic / Nonconservative)
 ---
 7. 빈출 템플릿: Rotating Hoop
 8. 빈출 템플릿: Rotating Triangle
 9. 빈출 템플릿: 2-DOF Chain System
 ---
 h. 도움말 (Help)
 0. 종료
---------------------------------------------
선택 >>
```

### 1.3 공통 입력 규칙

| 항목 | 규칙 |
|------|------|
| 행렬 | 행은 `;`, 열은 `,` 구분. 예: `9,0;0,9` → `[[9,0],[0,9]]` |
| 벡터 | `,` 구분. 예: `1,0` → `[1.0, 0.0]` |
| 선택 입력 | Enter를 누르면 기본값 사용 |
| 수식 | SymPy 문법. 분수는 `Rational(1,2)`, 거듭제곱은 `**` |

---

## 2. 메뉴 a — 자동 분류 (문제 유형 자동 판별)

**v2에서 새로 추가된 기능.** 문제 텍스트나 키워드를 입력하면 6가지 유형 중 가장 적합한 솔버를 자동으로 선택합니다.

### 사용 흐름

```
선택 >> a

=== 자동 문제 분류 ===
문제 텍스트 또는 키워드를 입력하세요.

입력: x''' + 3x'' + 6x' + 8x = f(t)

  분류 결과: [frf] FRF 분해 (전달함수)
  신뢰도: 50%
  매칭 키워드: ode

  [FRF 분해 (전달함수)] 메뉴를 실행하시겠습니까? (y/n) [y]:
```

### 지원하는 6가지 유형

| 유형 | 대표 키워드 / 패턴 | 연결 메뉴 |
|------|---------------------|:---------:|
| **tf** | true, false, T/F, 참, 거짓 | 1 |
| **frf** | transfer function, ODE, bode, G(s), 부분분수 | 2 |
| **lagrange** | Lagrange, T/V, 평형점, 선형화, virtual work | 3 |
| **modal** | eigenvalue, 고유값, mode shape, M/K 행렬 | 4 |
| **damping** | damping, Rayleigh, alpha/beta, zeta | 5 |
| **extended** | gyroscopic, flutter, divergence, 반대칭, bi-orthogonality | 6 |

### 입력 예시

```
x''' + 3x'' + 6x' + 8x = f(t)                    → frf
T = 1/2*m*rdot^2, V = -mgR*cos(theta), 평형점      → lagrange
M = [[9,0],[0,9]], K eigenvalue                     → modal
C = alpha*M + beta*K, damping ratio                 → damping
gyroscopic flutter divergence stability             → extended
이 문장이 참인지 거짓인지 판별                           → tf
```

---

## 3. 메뉴 1 — T/F 개념 판별 (v2)

**v2에서 대폭 개선.** 기존 키워드 매칭에 문장 구조 분석 + 충돌 검사가 추가되어, 시험형 함정 문제를 감지합니다.

### v1 vs v2 비교

| | v1 | v2 |
|--|:--:|:--:|
| 전체 정확도 | 50% (6/12) | **100% (12/12)** |
| FALSE 문항 정확도 | 14% (1/7) | **100% (7/7)** |
| 함정 감지 | 불가 | **부정어/조건/수치/범위 감지** |

### 사용 흐름

```
선택 >> 1

=== T/F 개념 판별 (v2) ===
T/F 문장을 입력하세요. 솔버가 구조 분석 후 판정합니다.

문장: Structural damping can be used for transient response analysis.
추가 키워드 (쉼표 구분, Enter to skip):
토픽 필터 (Enter to skip):
```

### 출력 예시

```
=======================================================
매칭 규칙: Structural damping (hysteretic damping) is valid for any arbitrary time-domain excitation.
규칙 verdict: FALSE
주의 포인트:
  - [trap_pattern] Structural damping 등가식은 harmonic excitation에서만 유효, transient/arbitrary에서는 사용 불가
추정 판정: FALSE
흔한 함정: Students apply structural damping to arbitrary time-domain problems; it is only valid under harmonic assumption.
유효 조건: harmonic excitation, steady-state harmonic
무효 조건: transient, arbitrary excitation, general time-domain, free vibration
-------------------------------------------------------
  ** 추정 판정: FALSE (충돌 1건 감지) **
=======================================================
```

### v2가 감지하는 함정 유형

| 유형 | 설명 | 예시 |
|------|------|------|
| **부정어 불일치** | 입력과 규칙의 부정어(not/cannot) 유무가 다름 | "cannot be diagonalized" vs 규칙 "can be diagonalized" |
| **함정 패턴** | 규칙에 정의된 특정 함정 정규식 매칭 | "both convolution" (주파수 영역은 multiplication) |
| **조건 충돌** | invalid_for에 해당하는 조건이 입력에 존재 | "free vibration" (convolution은 forced만) |
| **범위 충돌** | 도메인 키워드가 valid_for 범위 밖 | "transient" (structural damping은 harmonic만) |
| **수치 불일치** | 차원/수치가 규칙의 key_quantities와 다름 | "N-dimensional" (올바른 값: 2N) |
| **범위 과확장** | "both", "always" 등이 제한된 scope와 충돌 | "both free and forced" (forced만 유효) |

### 사용 팁

- 영어 문장을 그대로 입력하는 것이 가장 정확합니다
- 시험 문제의 T/F 문장을 복붙하세요
- 추가 키워드를 넣으면 관련 규칙 탐색 범위가 넓어집니다
- 토픽 필터로 특정 주제만 검색 가능: `damping`, `modal_analysis`, `gyroscopic_nonconservative` 등

### 사용 가능한 토픽 목록

```
damping, eigenvalue_problems, frf_convolution,
gyroscopic_nonconservative, linear_systems, modal_analysis,
natural_frequency, stability, state_equations,
structural_damping, virtual_work_lagrange
```

---

## 4. 메뉴 2 — FRF 분해 (전달함수)

### 개요

ODE를 입력하면 전달함수 G(s)를 구성하고, 부분분수 분해로 1차/2차 서브시스템의 파라미터를 추출합니다.

### 입력 흐름

```
선택 >> 2

=== FRF 분해 (전달함수) ===
  a) ODE 문자열  (예: x''' + 3x'' + 6x' + 8x = f(t))
  b) 계수 리스트 (예: 1,3,6,8)

입력 방식 (a/b) [a]: a
ODE 입력: x''' + 3x'' + 6x' + 8x = f(t)
분자 계수 (Enter for [1]):
강제 입력 f(t) (Enter to skip): 2cos(2)
```

### ODE 입력 방식

**방식 a) ODE 문자열:**

| 표기 | 예시 |
|------|------|
| 프라임(') | `x''' + 3x'' + 6x' + 8x = f(t)` |
| 캐럿 | `x^(3) + 3x^(2) + 6x^(1) + 8x = f(t)` |
| 소수 계수 | `x''' + 1.2x'' + 4.2x' + 4x = f(t)` |

**방식 b) 계수 리스트:**

```
계수 입력 (고차->저차, 쉼표 구분): 1,3,6,8
```

### 분자 계수

```
분자 계수 (Enter for [1]):
```

- 대부분 Enter (기본값 `[1]`)로 넘김
- 분자가 `s + 2` 같은 경우: `1,2` 입력

### 강제 입력 f(t)

Enter를 누르면 건너뛰고, 입력하면 해당 f(t)에 대한 응답을 계산합니다.

| # | 입력 | f(t) | 응답 계산 방법 |
|---|------|------|---------------|
| 1 | `impulse` 또는 `delta` | delta(t) | x(t) = L^{-1}{G(s)} |
| 2 | `3impulse` | 3*delta(t) | 스케일된 충격응답 |
| 3 | `step` 또는 `u(t)` | u(t) | x(t) = L^{-1}{G(s)/s} |
| 4 | `5step` | 5*u(t) | 스케일된 계단응답 |
| 5 | `ramp` 또는 `t` | t*u(t) | x(t) = L^{-1}{G(s)/s^2} |
| 6 | `exp(-2t)` | e^{-2t} | x(t) = L^{-1}{G(s)/(s+2)} |
| 7 | `cos` | cos(t) | x_ss = |G(j)|cos(t + angle G) |
| 8 | `3cos(2)` | 3cos(2t) | x_ss = 3|G(j2)|cos(2t + angle G) |
| 9 | `5sin(3)` | 5sin(3t) | x_ss = 5|G(j3)|sin(3t + angle G) |
| 10 | `2cos(1)+3sin(4)` | 중첩 | 중첩 원리로 합산 |

- `10` (숫자만) = `10step` (상수 = 스케일된 계단함수)
- impulse/step/ramp/exp: **완전한 시간응답** x(t) (역라플라스)
- cos/sin: **정상상태 응답** x_ss(t) (FRF 크기/위상)

### 출력 내용

| 항목 | 설명 |
|------|------|
| 전달함수 G(s) | 분모/분자 다항식 |
| 특성근 (poles) | 분모 다항식의 근 |
| 부분분수 분해 | 1차 + 2차 서브시스템 |
| 1차 서브시스템 | break frequency omega_b, 시정수 tau |
| 2차 서브시스템 | 고유진동수 omega_n, 감쇠비 zeta |
| G_i(jw) | 각 서브시스템의 jw 도메인 함수 (Re/Im 분리) |
| 특성값 테이블 | |G(0)|, |G(jw_b)|, |G(jw_n)|, |G(jinf)| |
| 강제응답 x(t) | f(t) 입력 시 응답 계산 |
| 안정성 | 모든 극점의 좌반면 위치 확인 |

### 출력 예시: jw 도메인

```
G1(jw) = [1] / [(12) + j*(6*omega)]
Re[G1(jw)] = 1/(3*(omega**2 + 4))
Im[G1(jw)] = -omega/(6*omega**2 + 24)
```

### 출력 예시: 특성값 테이블

```
G:
  |G(0)| = 0.125
  |G(jw_b)| = 0.176777  (angle -135.00 deg)
  |G(jw_n)| = 0.176777  (angle -135.00 deg)
  |G(jinf)| = 0
```

### 출력 예시: 강제응답

```
Forced response: f(t) = 2 cos(2t)
  |G(j2)| = 0.176777
  angle G(j2) = -135.00 deg
  x_ss(t) = 0.353553 cos(2t - 135.00 deg)
```

### 예제: 2025 기출

```
입력 방식 (a/b) [a]: a
ODE 입력: x''' + 3x'' + 6x' + 8x = f(t)
분자 계수 (Enter for [1]): (Enter)
강제 입력 f(t): 2cos(2)
```

---

## 5. 메뉴 3 — Lagrange 분석

### 개요

운동에너지 T와 위치에너지 V를 입력하면 Lagrange 방정식으로 EOM을 유도하고, 평형점 -> 선형화 -> 안정성 분석까지 수행합니다.

**v2에서 N-DOF 지원 추가.** DOF를 먼저 입력하여 1-DOF 또는 다자유도를 선택합니다.

### 사용 흐름

```
선택 >> 3

=== Lagrange 방정식 분석 ===
자유도 (DOF) [1]: 2

T (운동에너지) = Rational(1,2)*m1*q1dot**2 + Rational(1,2)*m2*q2dot**2
V (위치에너지) [0] = Rational(1,2)*k1*q1**2 + Rational(1,2)*k2*(q2-q1)**2
매개변수 이름: m1,m2,k1,k2
매개변수 수치값: m1=9,m2=9,k1=36,k2=36

좌표 1 이름 [q1]: q1
속도 1 이름 [q1dot]: q1dot
좌표 2 이름 [q2]: q2
속도 2 이름 [q2dot]: q2dot
비보존력 (예: q1=-c*q1dot,q2=0 / Enter to skip):
```

### 1-DOF 입력 (DOF=1)

```
자유도 (DOF) [1]: 1
T (운동에너지) = Rational(1,2)*m*R**2*(thetadot**2 + Omega**2*sin(theta)**2)
V (위치에너지) [0] = m*g*R*cos(theta)
일반화 좌표 이름 [q]: theta
일반화 속도 이름 [thetadot]: thetadot
매개변수 이름: m,g,R,Omega
매개변수 수치값: m=1,g=9.81,R=0.5,Omega=5
비보존력 Q_nc [0]: (Enter)
```

### N-DOF 입력 (DOF=2 이상)

```
자유도 (DOF) [1]: 2
T (운동에너지) = (한 줄에 모든 DOF의 T 합산)
V (위치에너지) = (한 줄에 모든 DOF의 V 합산)
매개변수 이름: m1,m2,k1,k2
매개변수 수치값: m1=9,m2=9,k1=36,k2=36
좌표 1 이름 [q1]: q1
속도 1 이름 [q1dot]: q1dot
좌표 2 이름 [q2]: q2
속도 2 이름 [q2dot]: q2dot
비보존력: q1=-c*q1dot,q2=0
```

### SymPy 입력 문법 규칙

| 표현 | 방법 | 예시 |
|------|------|------|
| 분수 1/2 | `Rational(1,2)` | `Rational(1,2)*m*v**2` |
| 거듭제곱 | `**` | `r**2`, `Omega**2` |
| 삼각함수 | `sin()`, `cos()` | `sin(theta)`, `cos(q1)` |
| 제곱근 | `sqrt()` | `sqrt(3)` |
| 원주율 | `pi` | `2*pi*r` |

> **주의**: `1/2`로 쓰면 Python 정수 나눗셈으로 `0`이 됩니다. 반드시 `Rational(1,2)` 사용!

### 출력 내용 (1-DOF)

| Step | 내용 |
|------|------|
| 1 | 일반화 좌표 정의 |
| 2-3 | T, V 파싱 |
| 4 | dT/d(qdot) |
| 5 | d/dt(dT/dqdot) (chain rule) |
| 6-7 | dT/dq, dV/dq |
| 8 | **EOM** 유도 |
| 9 | **평형점** (qdot=0, qddot=0) |
| 10-11 | **선형화**: m_eff, c_eff, k_eff 추출 |
| 12 | **안정성 판별** |

### 출력 내용 (N-DOF)

| Step | 내용 |
|------|------|
| 1 | N개 좌표 정의 |
| 2-3 | T, V 파싱 |
| 4 | 각 좌표의 **EOM** 유도 |
| 5 | **평형점** 계산 |
| 6 | **선형화**: M, C, K 행렬 자동 추출 |
| 7 | 수치 행렬 (parameter_values 입력 시) |
| 8 | **안정성**: K 고유값 기반 판별 |

### 예제: 2-DOF 직렬 스프링-질량

```
DOF: 2
T = Rational(1,2)*m1*q1dot**2 + Rational(1,2)*m2*q2dot**2
V = Rational(1,2)*k1*q1**2 + Rational(1,2)*k2*(q2-q1)**2
매개변수: m1,m2,k1,k2
수치값: m1=9,m2=9,k1=36,k2=36
```

**결과:**
- EOM_1: k1*q1 + k2*q1 - k2*q2 + m1*q1ddot = 0
- EOM_2: -k2*q1 + k2*q2 + m2*q2ddot = 0
- M = [[9,0],[0,9]], K = [[72,-36],[-36,36]]
- 안정성: STABLE (K 고유값 모두 양수)

### 예제: 1-DOF Rotating Hoop (2024 기출)

```
DOF: 1
T = Rational(1,2)*m*R**2*(thetadot**2 + Omega**2*sin(theta)**2)
V = m*g*R*cos(theta)
좌표: theta, 속도: thetadot
매개변수: m,g,R,Omega
수치값: m=1,g=9.81,R=0.5,Omega=5
```

### 주의사항

1. **좌표/속도 이름 일치** — T, V에서 쓴 이름과 입력하는 이름이 반드시 같아야 합니다
2. **매개변수 누락 금지** — T, V에 등장하는 기호를 모두 나열해야 합니다
3. **삼각함수 평형점** — 초월방정식은 SymPy가 모든 해를 찾지 못할 수 있습니다
4. **N-DOF 안정성** — 수치값(parameter_values)이 있어야 K 고유값 기반 안정성 판별이 가능합니다

---

## 6. 메뉴 4 — N-DOF Modal Analysis

### 개요

M(질량행렬)과 K(강성행렬)를 입력하면 고유값 문제 K u = lambda M u 를 풀어 고유진동수, 질량정규화 모드형상, 직교성 검증, 초기조건 응답을 계산합니다.

### 사용 흐름

```
선택 >> 4

=== N-DOF Modal Analysis ===
M (질량행렬): 9,0;0,9
K (강성행렬): 72,-36;-36,36
초기 변위 q(0): 1,0
초기 속도 qdot(0): 0,0
```

### 행렬 입력 형식

```
9,0;0,9        → [[9,0],[0,9]]
4,0;0,1        → [[4,0],[0,1]]
1,0,0;0,2,0;0,0,3  → 3x3 대각행렬
```

> M: 대칭 양정치 (positive-definite)
> K: 대칭 양반정치 (positive semi-definite)

### 출력 내용

| Step | 내용 |
|------|------|
| 1-2 | M, K 행렬 표시 |
| 3 | 고유값 lambda_i 계산 |
| 4 | **고유진동수** omega_i = sqrt(lambda_i) |
| 5 | **질량정규화 고유벡터** (u_i^T M u_i = 1) |
| 6 | 모달행렬 U = [u_1, u_2, ...] |
| 7 | **직교성 검증**: U^T M U = I |
| 8 | **직교성 검증**: U^T K U = Lambda |
| 9 | 모달 초기조건: eta(0) = U^T M q(0) |
| 10 | **모달 응답**: eta_r(t) = A cos(w_r t) + B sin(w_r t) |
| 11 | **물리 응답**: q(t) = U * eta(t) |

### M, K 행렬 조립 가이드

**직렬 체인: wall — k1 — [m1] — k2 — [m2] — k3 — wall**

```
M = [[m1, 0],        K = [[k1+k2, -k2  ],
     [0,  m2]]            [-k2,   k2+k3]]
```

**자유-자유: [m1] — k — [m2]**

```
M = [[m1, 0],        K = [[k,  -k],
     [0,  m2]]            [-k,  k]]
```

### 예제: 2025 기출

```
M: 9,0;0,9
K: 72,-36;-36,72
q(0): 1,1
qdot(0): 0,0
```

- lambda_1 = 4, lambda_2 = 12
- omega_1 = 2, omega_2 = 2*sqrt(3) rad/s
- u_1 = [1,1]/sqrt(18) (동위상), u_2 = [1,-1]/sqrt(18) (역위상)

---

## 7. 메뉴 5 — Proportional Damping Response

### 개요

M, K에 감쇠 정보를 추가하면 비례감쇠 여부를 판별하고 감쇠 모달 응답을 계산합니다.

### 감쇠 입력 3가지 방식

```
감쇠 입력 방식:
  a) C 행렬 직접 입력
  b) alpha, beta (Rayleigh damping: C = alpha*M + beta*K)
  c) 모달 감쇠비 zeta 직접 입력
```

**방식 a) C 행렬**
```
C (감쇠행렬): 2,-1;-1,2
```
- 자동으로 비례/비비례 판별
- 비비례이면 state-space 접근

**방식 b) alpha, beta**
```
alpha: 0.4
beta: 0.02
```
- C = alpha*M + beta*K 자동 구성
- 비례감쇠 보장
- 모달 감쇠비: zeta_r = (alpha + beta*omega_r^2) / (2*omega_r)

**방식 c) 모달 감쇠비**
```
모달 감쇠비: 0.05,0.1
```
- 각 모드의 zeta를 직접 지정

### 출력 내용

비례감쇠 시:
- 비감쇠 모달 분석 (고유진동수, 모드형상)
- 비례감쇠 판별: C = alpha*M + beta*K 확인
- **모달 감쇠비** zeta_r
- **감쇠 고유진동수** omega_d,r = omega_r * sqrt(1 - zeta_r^2)
- **감쇠 응답** q(t) = U * eta(t) (각 eta_r이 감쇠 진동)

비비례감쇠 시:
- State-space A 행렬 구성
- 복소 고유값에서 omega_n, zeta 추출

### alpha, beta 계산 가이드

두 모드의 감쇠비가 주어진 경우:

```
zeta_1 = (alpha + beta * omega_1^2) / (2 * omega_1)
zeta_2 = (alpha + beta * omega_2^2) / (2 * omega_2)
```

행렬로 풀면:
```
| 1/(2w1)   w1/2 | | alpha |   | zeta_1 |
| 1/(2w2)   w2/2 | | beta  | = | zeta_2 |
```

### 예제: 2025 기출 Rayleigh Damping

```
M: 9,0;0,9
K: 72,-36;-36,72
방식: b
alpha: 0.4
beta: 0.02
q(0): 1,1
qdot(0): 0,0
```

---

## 8. 메뉴 6 — 확장 (Gyroscopic / Nonconservative)

**v2에서 CLI 연결 완료.** M q'' + (C+G) q' + (K+H) q = 0 시스템을 해석합니다.

### 사용 흐름

```
선택 >> 6

=== 확장 시스템 (Gyroscopic / Nonconservative) ===
M q'' + (C+G) q' + (K+H) q = 0
행렬 입력: 행은 ';', 열은 ',' 구분

M (질량행렬): 1,0;0,1
K (강성행렬): 3,-1;-1,3
G (자이로 행렬, Enter to skip): 0,2;-2,0
H (비보존력 행렬, Enter to skip):
C (감쇠행렬, Enter to skip):
초기 변위 q(0) (Enter to skip):
초기 속도 qdot(0) (Enter to skip):
```

### 입력 행렬 정리

| 행렬 | 성질 | 필수 | 예시 |
|------|------|:----:|------|
| M | 대칭, 양정치 | O | `1,0;0,1` |
| K | 대칭 | O | `3,-1;-1,3` |
| G | **반대칭** (G^T = -G) | 선택 | `0,2;-2,0` |
| H | **반대칭** (H^T = -H) | 선택 | `0,1.5;-1.5,0` |
| C | 대칭 | 선택 | `0.5,0;0,0.5` |

> G와 H 중 최소 하나는 입력해야 이 메뉴의 의미가 있습니다.

### 출력 내용

| Step | 내용 |
|------|------|
| 1 | 행렬 성질 검증 (대칭/반대칭) |
| 2 | **시스템 분류**: pure_gyroscopic, damped_gyroscopic, pure_nonconservative 등 |
| 3 | State-space A 행렬 |
| 4 | **우고유벡터** (right): A phi = lambda phi |
| 5 | **좌고유벡터** (left): psi^T A = lambda psi^T |
| 6 | **쌍직교성 검증**: psi^T phi = I |
| 7 | **복소 모드**: omega_n, omega_d, zeta, sigma |
| 8 | **안정성 분류** |

### 안정성 분류

| 고유값 | 분류 |
|--------|------|
| 모두 Re(lambda) < 0 | **점근 안정** (asymptotically stable) |
| 모두 순허수, 단순근 | **자이로 안정** (gyroscopically stable) |
| Re > 0, Im != 0 | **Flutter 불안정** (진동 발산) |
| Re > 0, Im = 0 | **Divergence 불안정** (단조 발산) |

### 예제 1: 순수 자이로 시스템

```
M: 1,0;0,1
K: 5,0;0,5
G: 0,3;-3,0
```
- 시스템 유형: pure_gyroscopic
- 고유값: 순허수 -> 자이로스코픽 안정화

### 예제 2: 비보존력 (Flutter)

```
M: 1,0;0,1
K: 3,0;0,3
H: 0,1.5;-1.5,0
```
- H가 충분히 크면 flutter 불안정 발생

### 예제 3: 감쇠 + 자이로 복합

```
M: 2,0;0,2
K: 10,-2;-2,10
G: 0,4;-4,0
C: 0.5,0;0,0.5
```

### G, H 행렬 구성 규칙

자이로 행렬 G (반대칭):
```
G = [[  0,   g12, g13],
     [-g12,   0,  g23],
     [-g13, -g23,  0 ]]
```

비보존력 행렬 H (반대칭):
```
H = [[  0,   h12],
     [-h12,   0 ]]
```

---

## 9. 메뉴 7~9 — 빈출 템플릿

자주 출제되는 문제 유형을 미리 구성해 둔 템플릿입니다.
매개변수만 입력하면 전체 풀이가 자동으로 진행됩니다.

### 메뉴 7: Rotating Hoop

회전 후프(반지름 R, 각속도 Omega) 위의 비드(질량 m) 문제.

```
R (후프 반지름, m) [symbolic]: 0.5
m (질점 질량, kg) [symbolic]: 1
Omega (각속도, rad/s) [symbolic]: 5
g (중력가속도, m/s^2) [9.81]: 9.81
```

### 메뉴 8: Rotating Triangle

회전 삼각형 프레임 위의 질점 문제 (2025 기출).

```
m (질점 질량, kg) [symbolic]: 1
Omega (각속도, rad/s) [symbolic]: 5
g (중력가속도, m/s^2) [9.81]: 9.81
```

### 메뉴 9: 2-DOF Chain System

직렬 스프링-질량 시스템.

```
질량 [m1, m2]: 9,9
스프링 [k1, k2, k3]: 36,36,0
```

---

## 10. 결과 저장 및 그래프

### 저장

모든 메뉴 실행 후 결과 저장을 묻습니다:

```
결과를 파일로 저장하시겠습니까? (y/n): y
  -> Markdown: src/me551_solver/reports/20260413_123456_FRF_Decomposition.md
  -> Text:     src/me551_solver/reports/20260413_123456_FRF_Decomposition.txt
```

- Markdown: LaTeX 수식 블록 포함
- Text: 시험 답안 스타일 (Unicode)

### 그래프

FRF, Modal, Damping 메뉴에서는 그래프를 제공합니다:

```
그래프를 보시겠습니까? (y/n): y
그래프를 파일로 저장하시겠습니까? (y/n): y
```

- FRF: Bode plot (magnitude log-log + phase)
- Modal: 모드형상 시각화
- Damping: 시간응답 플롯

> matplotlib이 설치되어 있어야 합니다.

---

## 11. v2 변경사항 요약

### v1 -> v2 주요 변경

| 항목 | v1 | v2 |
|------|:--:|:--:|
| 버전 | 1.0 | **2.0** |
| **ConceptDB 정확도** | 50% | **100%** |
| **자동 분류** | 없음 | **메뉴 a (라우터)** |
| **Lagrange** | 1-DOF만 | **N-DOF 지원** |
| **Extended CLI** | 미연결 | **메뉴 6 연결** |
| Report 모드 | 1개 | **2개 (exam_style + exam_answer)** |
| 전체 테스트 | 133개 | **169개** |

### ConceptDB v2 상세

- 문장 구조 파서: 부정어/조건/수치/범위 감지
- 규칙별 trap_patterns: 정규식 기반 함정 패턴
- 규칙별 valid_for/invalid_for: 적용 범위 제한
- 충돌 검사: 6가지 유형의 구조적 차이 감지
- Stopword 필터: 기능어가 매칭을 왜곡하는 문제 해결

### 새로 추가된 파일

| 파일 | 설명 |
|------|------|
| `core/router.py` | 문제 유형 자동 분류기 |
| `tests/test_conceptdb_v2.py` | T/F 정확도 테스트 (14개) |
| `tests/test_router_e2e.py` | 라우터 + E2E 테스트 (22개) |

---

## 12. 트러블슈팅

### "Rational(1,2)" 관련

| 문제 | 원인 | 해결 |
|------|------|------|
| T에서 `1/2`를 쓰면 0이 됨 | Python 정수 나눗셈 | `Rational(1,2)` 사용 |
| 소수점 입력 시 거대 분수 | sympy.Rational(1.2) 문제 | v2에서 nsimplify로 수정됨 |

### 행렬 입력 관련

| 문제 | 원인 | 해결 |
|------|------|------|
| "행 수가 다릅니다" | `;` 누락 또는 `,` 개수 불일치 | 행렬 구조 확인 |
| M이 singular | 대각 원소에 0 | M은 양정치여야 함 |

### 좌표 이름 관련

| 문제 | 원인 | 해결 |
|------|------|------|
| "name 'rdot' is not defined" | 매개변수에 좌표를 넣음 | 좌표/속도는 매개변수에서 제외 |
| EOM이 0 = 0 | T, V에 쓴 이름과 좌표 이름이 다름 | 일치시킴 |

### 기타

| 문제 | 원인 | 해결 |
|------|------|------|
| matplotlib 관련 오류 | 미설치 | `pip install matplotlib` |
| 역라플라스 변환이 느림 | 복잡한 전달함수 | 정상 동작, 기다림 |
| T/F에서 "No rules matched" | 키워드 부족 | 추가 키워드 입력 또는 토픽 필터 사용 |
| Extended에서 G가 반대칭이 아님 | 입력 오류 | G[i][j] = -G[j][i] 확인 |
