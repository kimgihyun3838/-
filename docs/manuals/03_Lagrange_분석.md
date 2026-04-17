# 메뉴 3. Lagrange -> 평형점 -> 선형화 -> 안정성 — 사용 매뉴얼

## 개요

운동에너지 T와 위치에너지 V를 입력하면 Lagrange 방정식으로 EOM을 유도하고,
평형점 → 선형화 → 안정성 분석까지 자동 수행합니다. (1-DOF 전용)

---

## 실행 방법

```bash
cd C:/Users/김기현/Desktop/선형진동공학/src
python -m me551_solver
# 메뉴에서 3 입력
```

---

## 입력 흐름

### Step 1: 운동에너지 T

```
T (운동에너지) = Rational(1,2)*m*(rdot**2 + (Omega*r)**2)
```

**문법 규칙:**
- **SymPy 문법**으로 입력합니다
- 분수는 `Rational(분자, 분모)` 사용: `Rational(1,2)` = 1/2
- 거듭제곱은 `**` 사용: `r**2`, `Omega**2`
- 삼각함수: `sin(theta)`, `cos(theta)`, `tan(theta)`
- 좌표 속도는 이름 뒤에 `dot` 붙임: `rdot`, `thetadot`

**흔한 T 표현식:**

| 시스템 | T 표현식 |
|--------|----------|
| 병진 운동 | `Rational(1,2)*m*qdot**2` |
| 회전 운동 | `Rational(1,2)*I*thetadot**2` |
| 회전 후프 위 비드 | `Rational(1,2)*m*R**2*(thetadot**2 + Omega**2*sin(theta)**2)` |
| 회전 프레임 위 질점 | `Rational(1,2)*m*(rdot**2 + (Omega*r)**2)` |

### Step 2: 위치에너지 V

```
V (위치에너지) [0] = m*g*R*cos(theta)
```

- 기본값은 `0` (Enter 누르면 V=0)
- 중력 위치에너지 예: `m*g*R*cos(theta)`, `-m*g*r*sin(phi)`

**흔한 V 표현식:**

| 시스템 | V 표현식 |
|--------|----------|
| 없음 (수평 운동) | `0` |
| 단진자 | `-m*g*L*cos(theta)` |
| 수직 스프링 | `Rational(1,2)*k*q**2 + m*g*q` |
| 회전 후프 | `m*g*R*cos(theta)` |

### Step 3: 일반화 좌표 이름

```
일반화 좌표 이름 [q]: r
```

- 기본값은 `q`
- T, V에서 사용한 좌표 이름과 **일치**해야 함
- 예: `r`, `theta`, `q`, `x`

### Step 4: 일반화 속도 이름

```
일반화 속도 이름 [rdot]: rdot
```

- 기본값은 `(좌표이름)dot`
- T에서 사용한 속도 이름과 **일치**해야 함
- 예: `rdot`, `thetadot`, `qdot`, `xdot`

### Step 5: 매개변수 이름

```
매개변수 이름 (쉼표 구분, 예: m,g,R,Omega): m,g,R,Omega
```

- T, V에 등장하는 **모든** 기호 상수를 나열
- 좌표(`r`)와 속도(`rdot`)는 여기 적지 않음

### Step 6: 매개변수 수치값 (선택)

```
매개변수 수치값 (예: m=1,g=9.81 / Enter to skip): m=1,g=9.81,R=0.5,Omega=5
```

- Enter로 건너뛰면 **기호(symbolic) 해석**
- 수치를 넣으면 **수치 해석** (구체적인 숫자로 안정성 판별)
- 일부만 입력 가능: `g=9.81,R=1`

### Step 7: 비보존력 Q_nc

```
비보존력 Q_nc (Enter for 0) [0]: -c*rdot
```

- 기본값은 `0` (보존계)
- 점성 감쇠: `-c*rdot`, `-c*thetadot`
- 외력: `F`

---

## 출력 내용 (자동 계산 순서)

| Step | 내용 |
|------|------|
| 1 | 일반화 좌표 정의 |
| 2 | T 파싱 및 표시 |
| 3 | V 파싱 및 표시 |
| 4 | dT/d(qdot) 계산 |
| 5 | d/dt(dT/dqdot) 계산 (chain rule 적용) |
| 6 | dT/dq 계산 |
| 7 | dV/dq 계산 |
| 8 | **Lagrange 방정식 = EOM** 유도 |
| 9 | **평형점** 계산 (qdot=0, qddot=0으로 놓고 풀기) |
| 10 | 각 평형점에서 **섭동 도입** (q = q_eq + delta_q) |
| 11 | **선형화된 EOM**: m_eff * delta_qddot + c_eff * delta_qdot + k_eff * delta_q = 0 |
| 12 | **안정성 판별**: k_eff의 부호로 stable/unstable 분류 |

---

## 사용 예제

### 예제 1: Rotating Hoop (2024 기출)

회전 후프(반지름 R, 각속도 Omega) 위의 비드(질량 m) 문제.

```
T (운동에너지) = Rational(1,2)*m*R**2*(thetadot**2 + Omega**2*sin(theta)**2)
V (위치에너지) [0] = m*g*R*cos(theta)
일반화 좌표 이름 [q]: theta
일반화 속도 이름 [thetadot]: thetadot
매개변수 이름: m,g,R,Omega
매개변수 수치값: m=1,g=9.81,R=0.5,Omega=5
비보존력 Q_nc [0]: (Enter)
```

**결과:**
- EOM: m*R^2 * theta_ddot + m*g*R*sin(theta) - m*R^2*Omega^2*sin(theta)*cos(theta) = 0
- 평형점: theta = 0, theta = acos(g/(R*Omega^2))
- theta=0: Omega^2 > g/R일 때 불안정
- theta=acos(...): 안정

### 예제 2: Rotating Triangle (2025 기출)

회전 삼각형 프레임 위의 질점 문제.

```
T (운동에너지) = Rational(1,2)*m*(rdot**2 + (Omega*r)**2)
V (위치에너지) [0] = m*g*r
일반화 좌표 이름 [q]: r
일반화 속도 이름 [rdot]: rdot
매개변수 이름: m,g,Omega
매개변수 수치값: (Enter to skip)
비보존력 Q_nc [0]: (Enter)
```

### 예제 3: 단진자 (기본)

```
T (운동에너지) = Rational(1,2)*m*L**2*thetadot**2
V (위치에너지) [0] = -m*g*L*cos(theta)
일반화 좌표 이름 [q]: theta
일반화 속도 이름 [thetadot]: thetadot
매개변수 이름: m,g,L
매개변수 수치값: (Enter to skip)
비보존력 Q_nc [0]: (Enter)
```

---

## 주의사항

1. **1-DOF 전용** — 다자유도는 지원하지 않습니다
2. **좌표/속도 이름 일치** — T, V에서 쓴 이름과 Step 3, 4에서 입력하는 이름이 반드시 같아야 합니다
3. **Rational 사용** — `1/2`로 쓰면 Python 정수 나눗셈으로 `0`이 됩니다. 반드시 `Rational(1,2)` 사용
4. **매개변수 누락 금지** — T, V에 등장하는 기호를 모두 Step 5에 나열해야 정상 동작합니다
5. **삼각함수 평형점** — SymPy가 모든 해를 찾지 못할 수 있습니다 (특히 초월방정식)
