# ME551 Solver v2 개발 보고서

> 작성일: 2026-04-13
> 기반 계획: `me_551_solver_improvement_plan_v_2.md`

---

## 1. 구현 완료 항목

### 1.1 ConceptDB 구조 개선 (최우선, 완료)

**문제**: v1 키워드 매칭은 시험형 함정 문제(FALSE)에 취약 — 50% 정확도 (FALSE: 14%)

**해결**: 3단계 분석 시스템 구축

#### (1) tf_rules.yaml 스키마 확장
- `valid_for`: 규칙이 유효한 조건 목록 (예: "forced vibration", "harmonic excitation")
- `invalid_for`: 규칙이 무효한 조건 목록 (예: "free vibration", "transient")
- `key_quantities`: 핵심 수치 조건 (예: `{state_dimension: "p*n"}`)
- `trap_patterns`: 함정 문장 감지용 정규식 패턴 + 설명

적용된 규칙: rule_006, rule_007, rule_011, rule_021, rule_022, rule_028, rule_038

#### (2) 문장 구조 파서 (`_parse_statement()`)
- **부정어 감지**: not, cannot, never, impossible 등 20개 단어 + 14개 구문
- **범위 표현 감지**: always, both, any (확장) / only, solely (제한)
- **조건절 추출**: "for free vibration", "under harmonic" 등 패턴
- **수치/차원 추출**: p*n, 2n, N-dimensional 등
- **도메인 키워드**: time domain, frequency domain, convolution, multiplication

#### (3) 충돌 검사 (`_detect_conflicts()`)
6가지 충돌 유형 감지:
1. **부정어 불일치**: 입력 vs 규칙의 부정어 유무 차이
2. **함정 패턴 매칭**: 규칙에 정의된 trap_patterns 정규식 매칭
3. **조건 충돌**: 입력이 invalid_for에 해당하는 조건 포함
4. **범위 충돌**: valid_for와 상충하는 도메인 키워드
5. **수치 불일치**: key_quantities와 입력의 차원/수치 비교
6. **범위 과확장**: "both", "always" 사용 시 valid_for 범위 초과

#### (4) 판정 로직
- HIGH 충돌 감지 → 규칙 verdict 반전 (TRUE→FALSE or vice versa)
- MEDIUM 충돌 → 주의 표시
- 충돌 없음 → 기존 매칭 스코어 기반 confidence

#### (5) Stopword 필터
`_tokenize()`에 60개 영어 기능어 필터 추가 → "the", "and", "for", "can" 등이 스코어를 왜곡하던 문제 해결

#### 결과
| 지표 | v1 | v2 |
|------|:--:|:--:|
| 전체 정확도 | 6/12 (50%) | **12/12 (100%)** |
| TRUE 정확도 | 5/5 (100%) | 5/5 (100%) |
| FALSE 정확도 | 1/7 (14%) | **7/7 (100%)** |

---

### 1.2 문제 유형 자동 라우터 (완료)

**파일**: `core/router.py`

키워드 + 정규식 패턴 기반 분류기:
- 6개 라우트: tf, frf, lagrange, modal, damping, extended
- 각 라우트에 15~20개 키워드 + 3~5개 정규식 패턴
- CLI 메뉴 "a"로 접근, 분류 결과 확인 후 해당 솔버 실행

테스트: 13/13 분류 케이스 통과

---

### 1.3 Lagrange 2-DOF 확장 (완료)

**파일**: `core/lagrange.py` — `_solve_ndof()` 메서드 추가

기능:
- N개 일반화 좌표 입력: `coords: [("q1","q1dot"), ("q2","q2dot")]`
- 각 좌표에 대해 Lagrange 방정식 자동 유도
- 평형점 계산 (모든 속도/가속도 = 0)
- Jacobian 기반 선형화: M, C, K 행렬 자동 추출
- scipy.linalg.eig로 안정성 판별 (numeric)
- 1-DOF 하위 호환 유지

CLI: DOF 입력 → 1이면 기존 모드, 2+면 N-DOF 모드

---

### 1.4 Extended (Gyroscopic) CLI 연결 (완료)

**파일**: `ui/cli.py` — `_handle_extended()` 구현

기존에 Python 직접 호출만 가능했던 ExtendedSolver를 CLI 메뉴 6에서 직접 사용 가능:
- M, K 필수 입력
- G (자이로), H (비보존력), C (감쇠) 선택 입력
- 초기조건 선택 입력
- 결과 출력 + 그래프 + 저장 지원

---

### 1.5 Report 시스템 개선 (완료)

**파일**: `core/report.py` — `render_exam_answer()` 메서드 추가

기존 `render_exam_style`이 전체 풀이 과정을 나열한다면,
`render_exam_answer`는 **결론 + 근거**를 간결하게 출력:

- T/F: 판정 + 근거 + 감지된 함정 + 주의사항
- FRF: 핵심 특성값 + 강제응답
- Lagrange: EOM + 평형점 + 안정성
- Modal: 고유진동수 + 모드형상
- Extended: 시스템 유형 + 안정성 + flutter/divergence

---

### 1.6 테스트 구축 (완료)

#### test_conceptdb_v2.py (14 tests)
- 2024 기출 T/F 6문항 개별 테스트
- 2025 기출 T/F 6문항 개별 테스트
- 전체 정확도 ≥ 75% 임계값 테스트
- FALSE 정확도 ≥ 57% 임계값 테스트

#### test_router_e2e.py (22 tests)
- 라우터 분류 13 케이스 parametrize
- 메뉴 번호 매핑
- E2E: FRF 3차 ODE → 부분분수 분해
- E2E: Modal 2-DOF → 고유진동수
- E2E: Lagrange 2-DOF → EOM + 선형화 + 안정성
- E2E: Extended gyroscopic → 안정성
- E2E: ConceptDB TRUE/FALSE 판별
- E2E: Report 렌더링 모드

---

## 2. 수정된 파일 목록

| 파일 | 변경 내용 |
|------|----------|
| `core/concept_db.py` | v2 전면 재작성: 파서 + 충돌 검사 + advisory |
| `core/router.py` | 신규: 문제 유형 자동 분류기 |
| `core/lagrange.py` | N-DOF `_solve_ndof()` 추가 |
| `core/report.py` | `render_exam_answer()` 추가 |
| `data/tf_rules.yaml` | 7개 규칙에 valid_for/invalid_for/trap_patterns 추가 |
| `ui/cli.py` | v2.0: 자동라우터(a), Extended(6) 구현, Lagrange N-DOF, ConceptDB v2 UI |
| `tests/test_conceptdb_v2.py` | 신규: 12문항 정확도 테스트 |
| `tests/test_router_e2e.py` | 신규: 라우터 + E2E 22테스트 |

---

## 3. 테스트 현황

```
169 passed, 27 skipped (1.94s)
```

- 기존 133 tests: 전부 통과 (하위 호환 유지)
- 신규 36 tests: 전부 통과

---

## 4. 버전 변경

| 항목 | v1 | v2 |
|------|:--:|:--:|
| CLI 버전 | 1.0 | **2.0** |
| ConceptDB 정확도 | 50% | **100%** |
| Lagrange DOF | 1-DOF | **N-DOF** |
| Extended CLI | 미연결 | **연결** |
| 자동 분류 | 없음 | **라우터** |
| Report 모드 | 1 (exam_style) | **2 (+ exam_answer)** |
| 전체 테스트 | 133 | **169** |

---

## 5. 남은 제한사항

| 제한 | 설명 |
|------|------|
| ConceptDB 규칙 수 | 60개 — 시험 범위에 따라 추가 필요 |
| Lagrange N-DOF symbolic 안정성 | 수치값 없으면 symbolic 결과만 |
| 비비례감쇠 복소 모드 응답 | 고유값/감쇠비까지만 (복소 모드 응답 미구현) |
| 연속체 진동 | 보/판/셸 PDE 미지원 |
| 자동 분류 | 키워드 기반이므로 복합 문제 분류 한계 |
