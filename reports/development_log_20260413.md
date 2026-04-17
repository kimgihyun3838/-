# ME551 Solver 개발 진행 보고서

> 작성일: 2026-04-13

---

## 1. 프로젝트 개요

- **과목**: ME551 Linear Vibration Engineering (선형진동공학)
- **솔버 위치**: `src/me551_solver/`
- **실행**: `cd src && python -m me551_solver`
- **매뉴얼**: `docs/manuals/` 폴더

---

## 2. 솔버 전체 구조

### 2.1 코어 엔진 (core/)

| 모듈 | 기능 | 범용성 |
|------|------|:------:|
| `frf.py` | 전달함수 → 부분분수 분해, jw 도메인, 특성값, 강제응답 | O |
| `modal.py` | N-DOF K u = λ M u 고유값 문제, 모드형상, 시간응답 | O |
| `lagrange.py` | T, V → EOM 유도, 평형점, 선형화, 안정성 | 1-DOF만 |
| `damping.py` | 비례/비비례 감쇠 판별, 모달 감쇠 응답 | O |
| `stability.py` | 고유값 기반 안정성 분류 (stable/unstable/flutter/divergence) | O |
| `state_space.py` | 2차 → 1차 state-space 변환, 전이행렬, 강제응답 | O |
| `extended.py` | M q̈ + (C+G) q̇ + (K+H) q = 0 자이로/비보존력 해석 | O |
| `concept_db.py` | T/F 개념 판별 (60개 rule DB 기반 키워드 매칭) | O |
| `report.py` | Markdown / exam-style 텍스트 렌더링 | - |
| `visualizer.py` | Bode plot, pole-zero map, 모드형상, 시간응답 플롯 | - |

### 2.2 템플릿 (templates/)

| 모듈 | 문제 |
|------|------|
| `frf_decompose.py` | ODE 문자열/계수 → FRF 분해 편의 래퍼 |
| `rotating_hoop.py` | 회전 후프 위 비드 (2024 기출) |
| `rotating_triangle.py` | 회전 삼각형 프레임 (2025 기출) |
| `two_dof_chain.py` | 직렬 스프링-질량 2-DOF 체인 |

### 2.3 CLI 메뉴 (ui/cli.py)

```
1. T/F 개념 판별
2. FRF 분해 (전달함수)
3. Lagrange → 평형점 → 선형화 → 안정성
4. 2-DOF Modal Analysis
5. Proportional Damping Response
6. 확장 (Gyroscopic / Nonconservative) — CLI 미연결, Python 직접 호출
7-9. 빈출 템플릿 (Rotating Hoop / Triangle / 2-DOF Chain)
```

### 2.4 데이터

| 파일 | 내용 |
|------|------|
| `data/tf_rules.yaml` | T/F 개념 규칙 60개 (11개 토픽) |
| `data/lecture_index.yaml` | 강의 인덱스 |

### 2.5 테스트 (tests/)

| 파일 | 커버리지 |
|------|----------|
| `test_midterm_2024.py` | FRF(14), Rotating Hoop(11), 2-DOF Modal 등 38 tests |
| `test_midterm_2025.py` | FRF(8), Rotating Triangle(8), 2-DOF+Damping(12+) 등 31 tests |
| `test_report.py` | ReportEngine 렌더링 |
| `test_previous_midterms.py` | 기타 기출 |
| **전체** | **133 passed, 27 skipped** |

---

## 3. 금일 (2026-04-13) 수정 사항

### 3.1 소수 계수 파싱 버그 수정 (`core/frf.py`)

- **문제**: `sp.Rational(1.2)` → `5404319552844595/4503599627370496` (거대 분수)
  - 인수분해 실패 → 3차 다항식이 1차+2차로 분해되지 않음
- **수정**: `sp.nsimplify(c, rational=True)` 사용 → `1.2` → `6/5`
- **영향**: `_coeffs_to_poly()` 함수

### 3.2 G(jw) 도메인 출력 추가 (`core/frf.py`)

- 각 서브시스템(G₁, G₂)에 대해 s = jω 대입한 주파수 응답 함수 출력
- **실수부/허수부 분리**: Re[G(jω)], Im[G(jω)]
- **표준형 표시**: `[실수부 + j*허수부] / [실수부 + j*허수부]`
  - 기존 `sp.cancel()` 사용 시 `-I/(6*omega - 12*I)` 같은 비직관적 형태 → 제거

**출력 예시:**
```
G₁(jω) = [1] / [(12) + j*(6*omega)]
Re[G₁(jω)] = 1/(3*(omega**2 + 4))
Im[G₁(jω)] = -omega/(6*omega**2 + 24)
```

### 3.3 특성값 테이블 추가 (`core/frf.py`)

- G(전체), G₁, G₂ 각각에 대해 자동 계산:
  - **|G(0)|** — DC gain
  - **|G(jω_b)|** — break frequency에서의 크기와 위상
  - **|G(jω_n)|** — natural frequency에서의 크기와 위상
  - **|G(j∞)|** — 고주파 극한

**출력 예시:**
```
G:
  |G(0)| = 0.125
  |G(jω_b)| = 0.176777  (∠ -135.00°)
  |G(jω_n)| = 0.176777  (∠ -135.00°)
  |G(j∞)| = 0

G1:
  |G1(0)| = 0.0833333
  |G1(jω_b)| = 0.0589256  (∠ -45.00°)
  ...
```

### 3.4 강제응답 x(t) 계산 추가 (`core/frf.py` + `ui/cli.py`)

#### 지원 입력 형식

| 입력 | f(t) | 계산 방법 |
|------|------|----------|
| `impulse`, `delta` | δ(t) | x(t) = L⁻¹{G(s)} |
| `3impulse` | 3δ(t) | 스케일된 충격응답 |
| `step`, `u(t)` | u(t) | x(t) = L⁻¹{G(s)/s} |
| `5step` | 5u(t) | 스케일된 계단응답 |
| `ramp`, `t` | t·u(t) | x(t) = L⁻¹{G(s)/s²} |
| `2t` | 2t·u(t) | 스케일된 램프응답 |
| `exp(-2t)` | e⁻²ᵗ | x(t) = L⁻¹{G(s)/(s+2)} |
| `cos`, `3cos(2)` | cos(t), 3cos(2t) | x_ss = \|G(jω₀)\|F₀cos(ω₀t + ∠G) |
| `5sin(3)` | 5sin(3t) | 정상상태 응답 |
| `2cos(1)+3sin(4)` | 중첩 | 중첩 원리 |
| `10` (숫자만) | 10·u(t) | = `10step` |

#### 동작 원리

- **impulse/step/ramp/exp**: SymPy `inverse_laplace_transform`으로 **완전한 시간응답 x(t)** 계산
- **cos/sin 조화입력**: 주파수 응답 함수로 **정상상태 응답 x_ss(t)** 계산
  - 각 서브시스템(G₁, G₂)의 기여도(크기/위상)도 함께 표시

**출력 예시 (step):**
```
f(t) = 1.0 u(t)
X(s) = G(s)/s = 1/(s⁴ + 3s³ + 6s² + 8s)
x(t) = 1/8 - exp(-2t)/12 - √15·exp(-t/2)·sin(√15·t/2)/40
       - exp(-t/2)·cos(√15·t/2)/24
```

**출력 예시 (harmonic):**
```
f(t) = 2 cos(2t)
|G(j2)| = 0.176777,  ∠G(j2) = -135.00°
x_ss(t) = 0.353553 cos(2t - 135.00°)
```

### 3.5 Bode plot 오류 수정 (`templates/frf_decompose.py`)

- **문제**: ODE 문자열 입력 시 `given`에 `coefficients` 키가 없어 visualizer가 Bode 데이터를 계산 못 함
- **수정**: `given["coefficients"] = parsed_coeffs` 추가

### 3.6 Bode plot y축 변경 (`core/visualizer.py`)

- **변경**: y축 dB 스케일 → 선형 magnitude (log-log)
- y축 라벨: `Magnitude (dB)` → `Magnitude |G(jω)|`
- `semilogx` → `loglog`
- 서브시스템 기여도 플롯도 동일하게 변경

### 3.7 사용 매뉴얼 작성 (`docs/manuals/`)

| 파일 | 내용 |
|------|------|
| `02_FRF_분해.md` | 입력법, jω 출력, 특성값 테이블, 강제응답 (10가지 예시), 예제 |
| `03_Lagrange_분석.md` | T/V 입력 문법, SymPy 규칙, Rational 주의사항, 예제 3개 |
| `04_Modal_Analysis.md` | M/K 행렬 입력, N-DOF 지원, 조립 가이드, 예제 4개 |
| `05_Proportional_Damping.md` | 감쇠 입력 3방식, 비례/비비례 판별, α-β 계산 가이드 |
| `06_Extended_Gyroscopic.md` | Python 직접 호출법, G/H 행렬 구성, flutter/divergence |

---

## 4. 수정된 파일 목록

| 파일 | 변경 내용 |
|------|----------|
| `src/me551_solver/core/frf.py` | nsimplify 수정, jω 출력, 특성값 테이블, 강제응답 |
| `src/me551_solver/core/visualizer.py` | Bode y축 dB→magnitude, loglog |
| `src/me551_solver/templates/frf_decompose.py` | coefficients 전달, forcing 전달 |
| `src/me551_solver/ui/cli.py` | forcing 파서, CLI 안내 문구 |
| `docs/manuals/02_FRF_분해.md` | 신규 작성 + 3회 업데이트 |
| `docs/manuals/03_Lagrange_분석.md` | 신규 작성 |
| `docs/manuals/04_Modal_Analysis.md` | 신규 작성 |
| `docs/manuals/05_Proportional_Damping.md` | 신규 작성 |
| `docs/manuals/06_Extended_Gyroscopic.md` | 신규 작성 |

---

## 5. 테스트 현황

```
133 passed, 27 skipped (0.15s)
```

기존 테스트 전체 통과 확인. 신규 기능은 수동 검증 완료.

---

## 6. 알려진 제한사항

| 제한 | 설명 |
|------|------|
| Lagrange 1-DOF만 | 다자유도 라그랑주는 수동 유도 필요 |
| Extended (메뉴 6) CLI 미연결 | Python 직접 호출로 사용 가능 |
| 비선형 해석 불가 | 선형화까지만 지원 |
| 연속체 진동 없음 | 보/판/셸 PDE 미지원 |
| 시변 계수 미지원 | Mathieu 방정식 등 불가 |
| 비비례감쇠 | 고유값/감쇠비 추출까지만 (복소 모드 응답 미구현) |
| 역라플라스 | 복잡한 전달함수에서 SymPy 변환 시간이 길 수 있음 |

---

## 7. ConceptDB T/F 정확도 테스트 결과

> 상세 리포트: `reports/ConceptDB_TF_Accuracy_Report.md`

### 테스트 대상
- 2024 기출 T/F 6문항 + 2025 기출 T/F 6문항 = **총 12문항**
- (2020-2023 기출에는 T/F 문제 없음)

### 결과: **6/12 (50%)**

| 연도 | 정답 | 정확도 |
|------|------|--------|
| 2024 | 4/6 | 67% |
| 2025 | 2/6 | 33% |
| **전체** | **6/12** | **50%** |

| 문항 유형 | 정답 | 정확도 |
|-----------|------|--------|
| TRUE 문항 | 5/5 | 100% |
| FALSE 문항 | 1/7 | **14%** |

### 오답 원인: 3가지 패턴

1. **참인 개념에 거짓 조건 추가** (3건) — "free vibration에도", "transient에도", "convolution" vs "multiplication"
2. **수치/차원 변조** (2건) — "p*n"을 "p" 또는 "N"으로 바꿔치기
3. **부정어 무시** (1건) — "cannot be diagonalized"의 "cannot" 무시

### 근본 한계

키워드 매칭은 **"가장 비슷한 rule 찾기"**만 수행하고, **"입력 문장이 그 rule과 같은 의미인지"**는 판단하지 못함.
시험 T/F는 정확히 이 약점을 공격하는 구조 (참 개념 + 미세 변형 → 거짓).

---

## 8. 향후 작업 후보

- [ ] 메뉴 6 (Extended) CLI 연결
- [ ] Lagrange 다자유도 확장
- [ ] 비비례감쇠 복소 모드 응답
- [ ] FRF Bode plot에 특성값 점 표시 (|G(0)|, |G(jω_n)| 등)
- [ ] 강제응답 시간 플롯 자동 생성
- [ ] ConceptDB 개선: 부정어 감지, 수치 비교, 조건절 파싱 등
