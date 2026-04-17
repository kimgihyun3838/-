# Project Template — 폴더 구조 INDEX

> 이 파일은 프로젝트 전체 폴더 구조를 설명합니다.
> Claude에게 "이 INDEX.md를 읽고 프로젝트 구조를 파악해줘"라고 하면 됩니다.

## 폴더 구조

```
my-project/
│
├── CLAUDE.md              ★ AI 업무 매뉴얼 (Claude가 자동으로 읽음)
├── plan.md                ★ 프로젝트 전체 계획 (Phase별, Gate 포함)
├── task.md                ★ 현재 작업 체크리스트 (세션별)
├── progress.md            ★ 진행 상황 스냅샷 (/clear 전 필수 업데이트)
├── lessons.md               실수에서 배운 교훈 모음
├── decisions.md              설계 결정 기록 (ADR-lite)
├── INDEX.md                  이 파일 — 폴더 구조 설명
│
├── src/                      소스 코드
│   ├── app.py                  메인 앱
│   ├── data_loader.py          데이터 로딩 모듈
│   └── visualizer.py           그래프 생성 모듈
│
├── tests/                    테스트 코드
│   ├── test_data_loader.py
│   └── test_visualizer.py
│
├── data/                     데이터 파일
│   ├── sample.csv              샘플 데이터
│   └── README.md               데이터 설명
│
├── references/               ★ 외부 조사 자료 (Claude Code 밖에서 수집)
│   ├── INDEX.md                참고 자료 목차
│   ├── streamlit_guide.md      웹에서 조사한 Streamlit 가이드
│   ├── plotly_examples.md      Plotly 사용법 정리
│   └── related_paper.pdf       관련 논문 PDF
│
├── papers/                   논문 관련 (연구 프로젝트인 경우)
│   ├── PLAN.md                 논문 집필 계획
│   ├── paper_main.tex          논문 본문
│   ├── references.bib          참고문헌 (BibTeX)
│   └── figures/                논문 그림
│       ├── fig01_overview.pdf
│       └── fig_src/            그림 생성 스크립트
│
├── research/                 연구 탐색 (실험/분석)
│   ├── INDEX.md                연구 주제 목록 + 상태
│   └── R1_data_exploration/    연구 주제별 폴더
│       ├── PLAN.md               이 연구의 계획
│       ├── analysis.py           분석 코드
│       ├── data/                 실험 데이터
│       └── figures/              결과 그래프
│
├── reports/                  보고서 출력물
│   ├── report_v1.html
│   └── report_v1.pdf
│
├── docs/                     문서 아카이브
│   ├── tasks/                  세션별 task 파일 (영구 보존)
│   │   ├── task_S001.md
│   │   └── task_S002.md
│   └── SESSION_HISTORY.md      세션 기록 (append-only)
│
├── requirements.txt          Python 패키지 목록
└── .gitignore                Git 제외 파일 설정
```

## 파일 역할 요약

### ★ 핵심 문서 (항상 최신 상태 유지)
| 파일 | 역할 | 언제 업데이트 |
|------|------|-------------|
| `CLAUDE.md` | AI 규칙/금지사항/워크플로우 | 새 규칙 발견 시 |
| `plan.md` | 전체 계획 (Phase + Gate) | 계획 변경 시 |
| `task.md` | 현재 세션 체크리스트 | 작업 완료 시 |
| `progress.md` | 진행 상황 스냅샷 | **`/clear` 전 필수** |

### 보조 문서
| 파일 | 역할 | 언제 업데이트 |
|------|------|-------------|
| `lessons.md` | 실수 기록 → 재발 방지 | 문제 발생 시 |
| `decisions.md` | 설계 결정 + 대안 + 사유 | 결정 시 |
| `references/INDEX.md` | 참고 자료 목록 | 자료 추가 시 |
| `research/INDEX.md` | 연구 주제 목록 + 상태 | 연구 추가/완료 시 |
| `docs/SESSION_HISTORY.md` | 세션별 상세 기록 | 세션 종료 시 |

### 폴더 용도
| 폴더 | 용도 | 비고 |
|------|------|------|
| `src/` | 소스 코드 | 메인 구현 |
| `tests/` | 테스트 코드 | pytest |
| `data/` | 데이터 파일 | CSV, JSON 등 |
| `references/` | 외부 조사 자료 | **Claude Code 밖에서 수집** |
| `papers/` | 논문 관련 | tex, bib, figures |
| `research/` | 연구 탐색 | R1, R2... 주제별 |
| `reports/` | 보고서 출력물 | HTML, PDF, DOCX |
| `docs/tasks/` | 세션별 task 아카이브 | 영구 보존 |
