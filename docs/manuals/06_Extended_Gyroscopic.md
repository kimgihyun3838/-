# 메뉴 6. 확장 (Gyroscopic / Nonconservative) — 사용 매뉴얼

## 개요

일반 운동방정식 M q_ddot + (C + G) q_dot + (K + H) q = 0 을 해석합니다.

- **G**: 자이로스코픽 행렬 (반대칭, G^T = -G)
- **H**: 비보존력/순환력 행렬 (반대칭, H^T = -H)

좌/우 고유벡터, 쌍직교성(bi-orthogonality) 검증, 복소 모드 분석,
flutter/divergence 불안정성 판별을 수행합니다.

---

## 실행 방법

```bash
cd C:/Users/김기현/Desktop/선형진동공학/src
python -m me551_solver
# 메뉴에서 6 입력
```

> **참고:** CLI 메뉴에서는 "v2에서 지원 예정"이라고 표시되지만,
> `ExtendedSolver` 엔진 자체는 이미 구현되어 있습니다.
> 프로그래밍 방식으로는 사용 가능합니다 (아래 "프로그래밍 사용법" 참조).

---

## 프로그래밍 사용법 (Python 코드에서 직접 호출)

CLI 메뉴가 아직 연결되지 않았으므로, Python에서 직접 호출합니다:

```python
import numpy as np
from me551_solver.core.extended import ExtendedSolver
from me551_solver.core.report import ReportEngine

solver = ExtendedSolver()

params = {
    "M": [[4, 0], [0, 4]],          # 질량행렬 (대칭 양정치)
    "K": [[8, 0], [0, 8]],          # 강성행렬 (대칭)
    "G": [[0, 2], [-2, 0]],         # 자이로스코픽 (반대칭)
    # "H": [[0, 1], [-1, 0]],       # 비보존력 (반대칭, 선택)
    # "C": [[0.1, 0], [0, 0.1]],    # 감쇠 (대칭, 선택)
}

result = solver.solve(params)

# 결과 출력
print(ReportEngine.render_exam_style(result))

# 파일 저장
ReportEngine.save_markdown(result, "gyroscopic_result.md")
```

실행:
```bash
cd C:/Users/김기현/Desktop/선형진동공학/src
python -c "위 코드"
```

또는 스크립트 파일로 저장 후 실행.

---

## 입력 파라미터

| 파라미터 | 타입 | 필수 | 설명 |
|----------|------|:----:|------|
| `M` | 2D list / ndarray | O | 질량행렬 (n x n, 대칭 양정치) |
| `K` | 2D list / ndarray | O | 강성행렬 (n x n, 대칭) |
| `G` | 2D list / ndarray | 선택 | 자이로스코픽 행렬 (n x n, 반대칭: G^T = -G) |
| `H` | 2D list / ndarray | 선택 | 비보존력 행렬 (n x n, 반대칭: H^T = -H) |
| `C` | 2D list / ndarray | 선택 | 감쇠행렬 (n x n, 대칭) |

> G와 H 중 최소 하나는 입력해야 의미가 있습니다.
> 둘 다 없으면 일반 모달 분석(메뉴 4)과 동일합니다.

### 행렬 성질 요약

| 행렬 | 성질 | 예시 (2x2) |
|------|------|------------|
| M (질량) | 대칭, 양정치 | `[[4,0],[0,4]]` |
| K (강성) | 대칭 | `[[8,-2],[-2,8]]` |
| G (자이로) | 반대칭 (G^T = -G) | `[[0,g],[-g,0]]` |
| H (비보존) | 반대칭 (H^T = -H) | `[[0,h],[-h,0]]` |
| C (감쇠) | 대칭 | `[[c1,0],[0,c2]]` |

---

## 출력 내용

| Step | 내용 |
|------|------|
| 1 | 입력 행렬 표시 및 대칭/반대칭 성질 검증 |
| 2 | **시스템 분류**: pure gyroscopic, damped gyroscopic, nonconservative 등 |
| 3 | State-space 정식화: A 행렬 구성 |
| 4 | **우고유벡터** (right eigenvectors): A phi_r = lambda_r phi_r |
| 5 | **좌고유벡터** (left eigenvectors): psi^T A = lambda psi^T |
| 6 | **쌍직교성 검증**: psi_i^T phi_j 관계 확인 |
| 7 | **복소 모드형상**: 각 모드의 진동 패턴 |
| 8 | **안정성 분류**: 고유값 위치 기반 판별 |

### 안정성 분류 기준

| 고유값 위치 | 분류 |
|-------------|------|
| 모두 Re(lambda) < 0 | 점근 안정 (asymptotically stable) |
| 모두 Re(lambda) <= 0, 순허수 단순근 | 한계 안정 (marginally stable, 자이로 안정화) |
| Re(lambda) > 0, Im(lambda) != 0 | **Flutter 불안정** (진동 발산) |
| Re(lambda) > 0, Im(lambda) = 0 | **Divergence 불안정** (단조 발산) |

---

## 사용 예제

### 예제 1: 순수 자이로스코픽 시스템

```python
params = {
    "M": [[1, 0], [0, 1]],
    "K": [[5, 0], [0, 5]],
    "G": [[0, 3], [-3, 0]],    # 자이로 효과
}
result = solver.solve(params)
```

- G가 K의 불안정성을 보상하여 **자이로스코픽 안정화** 가능

### 예제 2: 비보존력 시스템 (Flutter)

```python
params = {
    "M": [[1, 0], [0, 1]],
    "K": [[3, 0], [0, 3]],
    "H": [[0, 1.5], [-1.5, 0]],    # 순환력
}
result = solver.solve(params)
```

- H가 충분히 크면 **flutter 불안정** 발생

### 예제 3: 감쇠 + 자이로 복합 시스템

```python
params = {
    "M": [[2, 0], [0, 2]],
    "K": [[10, -2], [-2, 10]],
    "G": [[0, 4], [-4, 0]],
    "C": [[0.5, 0], [0, 0.5]],
}
result = solver.solve(params)
```

### 예제 4: 3-DOF 자이로스코픽

```python
params = {
    "M": [[1,0,0],[0,1,0],[0,0,1]],
    "K": [[4,-1,0],[-1,4,-1],[0,-1,4]],
    "G": [[0,2,0],[-2,0,1],[0,-1,0]],
}
result = solver.solve(params)
```

---

## 자이로스코픽 행렬 G 구성 가이드

자이로스코픽 효과가 있는 시스템의 G는 반대칭이어야 합니다:

```
G = [[  0,   g12, g13],
     [-g12,   0,  g23],
     [-g13, -g23,  0 ]]
```

**물리적 의미:**
- 회전체(로터, 자이로스코프)에서 코리올리 효과
- G가 하는 일(work) = 0 (에너지 보존)
- G는 고유진동수를 분리(split)시킴

---

## 비보존력 행렬 H 구성 가이드

순환력(circulatory force)에서 발생하는 H도 반대칭:

```
H = [[  0,   h12],
     [-h12,   0 ]]
```

**물리적 의미:**
- 종동력(follower force) 등 방향이 변하는 하중
- H가 충분히 크면 flutter 불안정 유발
- 감쇠가 오히려 불안정을 촉진할 수 있음 (destabilizing effect of damping)

---

## 주의사항

1. **CLI 메뉴 미연결** — 현재 CLI에서 6번을 선택하면 "v2 예정" 메시지만 출력됩니다. Python 코드로 직접 호출해야 합니다.
2. **반대칭 검증** — G와 H가 반대칭이 아니면 경고가 표시됩니다
3. **M 양정치 필수** — M의 역행렬이 필요하므로 특이행렬이면 오류 발생
4. **복소 고유값** — 이 솔버의 고유값은 일반적으로 복소수입니다 (대칭 고유값 문제가 아님)
5. **쌍직교성** — 좌/우 고유벡터 간의 직교 관계로, 대칭 시스템의 U^T M U = I 의 일반화
