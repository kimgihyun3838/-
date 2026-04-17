# 메뉴 5. Proportional Damping Response — 사용 매뉴얼

## 개요

M, K 행렬에 감쇠 정보(C 행렬 또는 alpha/beta 또는 모달 감쇠비)를 추가로 입력하면,
비례감쇠(C = alpha*M + beta*K) 여부를 확인하고, 감쇠 모달 응답을 계산합니다.
비비례감쇠인 경우에도 state-space 접근으로 고유진동수/감쇠비를 추출합니다.

---

## 실행 방법

```bash
cd C:/Users/김기현/Desktop/선형진동공학/src
python -m me551_solver
# 메뉴에서 5 입력
```

---

## 입력 흐름

### Step 1: 질량행렬 M

```
M (질량행렬, 예: 9,0;0,9): 9,0;0,9
```

행은 `;`, 열은 `,`으로 구분.

### Step 2: 강성행렬 K

```
K (강성행렬, 예: 72,-36;-36,72): 72,-36;-36,72
```

### Step 3: 감쇠 입력 방식 선택

```
감쇠 입력 방식:
  a) C 행렬 직접 입력
  b) alpha, beta (Rayleigh damping: C = alpha*M + beta*K)
  c) 모달 감쇠비 zeta 직접 입력
감쇠 입력 방식 (a/b/c) [b]:
```

---

### 방식 a) C 행렬 직접 입력

```
C (감쇠행렬): 2,-1;-1,2
```

- 솔버가 자동으로 **비례감쇠 여부를 판별**합니다
  - Least-squares: C = alpha*M + beta*K 에서 alpha, beta 추정
  - Commutativity check: M^{-1}C * M^{-1}K = M^{-1}K * M^{-1}C ?
- 비례감쇠이면: 모달 분해 → 감쇠 응답
- 비비례감쇠이면: state-space 고유값 분석

### 방식 b) alpha, beta 입력 (Rayleigh Damping)

```
alpha: 0.5
beta: 0.02
```

- C = alpha*M + beta*K 로 감쇠행렬 자동 구성
- 비례감쇠가 **보장**됨
- 모달 감쇠비: zeta_r = (alpha + beta * omega_r^2) / (2 * omega_r)

### 방식 c) 모달 감쇠비 직접 입력

```
모달 감쇠비 (쉼표 구분, 예: 0.05,0.1): 0.05,0.1
```

- 각 모드의 감쇠비를 직접 지정
- 모드 수만큼 입력 (2-DOF면 2개, 3-DOF면 3개)

### Step 4: 초기 변위 (선택)

```
초기 변위 q(0) (Enter to skip): 1,0
```

### Step 5: 초기 속도 (선택)

```
초기 속도 qdot(0) (Enter to skip): 0,0
```

---

## 출력 내용

### 비례감쇠인 경우

| Step | 내용 |
|------|------|
| 1 | 비감쇠 모달 분석 (ModalSolver 호출) — 고유진동수, 모드형상 |
| 2 | C 행렬 구성 또는 표시 |
| 3 | **비례감쇠 판별**: C = alpha*M + beta*K 확인, 잔차/교환성 검사 |
| 4 | **모달 감쇠비**: zeta_r = (alpha + beta*omega_r^2) / (2*omega_r) |
| 5 | **감쇠 고유진동수**: omega_d,r = omega_r * sqrt(1 - zeta_r^2) |
| 6 | 모달 방정식: eta_ddot_r + 2*zeta_r*omega_r*eta_dot_r + omega_r^2*eta_r = 0 |
| 7 | **모달 감쇠 응답**: eta_r(t) = e^(-sigma_r*t) [A_r cos(omega_d,r*t) + B_r sin(omega_d,r*t)] |
| 8 | **물리 감쇠 응답**: q(t) = U * eta(t) |

### 비비례감쇠인 경우

| Step | 내용 |
|------|------|
| 1 | 비감쇠 모달 분석 |
| 2 | C 행렬 표시 |
| 3 | 비례감쇠 판별: **NOT PROPORTIONAL** |
| 4 | State-space 고유값 분석: A = [[0,I],[-M^{-1}K, -M^{-1}C]] |
| 5 | 복소 고유값에서 omega_n, zeta, omega_d 추출 |

---

## 사용 예제

### 예제 1: 2025 기출 — Rayleigh Damping

```
M (질량행렬): 9,0;0,9
K (강성행렬): 72,-36;-36,72
감쇠 입력 방식 (a/b/c) [b]: b
alpha: 0.4
beta: 0.02
초기 변위 q(0): 1,1
초기 속도 qdot(0): 0,0
```

**결과:**
- omega_1 = 2, omega_2 = 2*sqrt(3) rad/s
- zeta_1 = (0.4 + 0.02*4) / (2*2) = 0.12
- zeta_2 = (0.4 + 0.02*12) / (2*2*sqrt(3)) = ...
- omega_d,1 = 2*sqrt(1 - 0.12^2)
- 각 q_i(t)에 대한 감쇠 응답식 출력

### 예제 2: C 행렬 직접 — 비례감쇠 자동 판별

```
M (질량행렬): 2,0;0,1
K (강성행렬): 6,-2;-2,2
감쇠 입력 방식 (a/b/c) [b]: a
C (감쇠행렬): 1.4,-0.4;-0.4,0.6
초기 변위 q(0): (Enter)
초기 속도 qdot(0): (Enter)
```

솔버가 C = alpha*M + beta*K 여부를 자동 검사합니다.

### 예제 3: 모달 감쇠비 직접 지정

```
M (질량행렬): 9,0;0,9
K (강성행렬): 72,-36;-36,72
감쇠 입력 방식 (a/b/c) [b]: c
모달 감쇠비: 0.05,0.1
초기 변위 q(0): 1,0
초기 속도 qdot(0): 0,0
```

---

## alpha, beta 계산 가이드

시험에서 alpha, beta가 주어지지 않고 두 모드의 감쇠비가 주어지는 경우:

```
zeta_1 = (alpha + beta * omega_1^2) / (2 * omega_1)
zeta_2 = (alpha + beta * omega_2^2) / (2 * omega_2)
```

이 연립방정식을 풀면:
```
| 1/(2*omega_1)   omega_1/2 | | alpha |   | zeta_1 |
| 1/(2*omega_2)   omega_2/2 | | beta  | = | zeta_2 |
```

---

## 주의사항

1. **감쇠 입력 3가지 중 하나만** 선택 — C 행렬, alpha+beta, zeta 중 1개
2. **비비례감쇠 시** 복소 모드형상이 필요하지만 현재 버전에서는 고유값/감쇠비 추출까지만 지원
3. **과감쇠(zeta >= 1)** 모드는 별도 처리됨 (omega_d = 0, overdamped 표시)
4. **C 행렬 방식 선택 시** 자동으로 비례/비비례 판별 → 적절한 해석 방법 선택
