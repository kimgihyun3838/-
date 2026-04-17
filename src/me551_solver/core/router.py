"""문제 유형 자동 라우터 — 입력 텍스트에서 문제 유형을 분류하여 적절한 엔진 호출.

키워드 기반 분류기로, 입력 텍스트를 분석하여 어떤 솔버 모듈을 사용할지 결정한다.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any


@dataclass
class RouteResult:
    """라우터 분류 결과."""
    route: str                     # 라우트 ID (e.g., "frf", "lagrange", ...)
    label: str                     # 사람 읽기용 라벨
    confidence: float              # 0.0 ~ 1.0
    matched_keywords: list[str] = field(default_factory=list)
    suggestions: list[str] = field(default_factory=list)  # 대안 라우트


# ---------------------------------------------------------------------------
# 라우트 정의: (route_id, label, keywords, weight_boost_patterns)
# ---------------------------------------------------------------------------
_ROUTES: list[dict[str, Any]] = [
    {
        "id": "tf",
        "label": "T/F 개념 판별",
        "menu": "1",
        "keywords": [
            "true", "false", "t/f", "참", "거짓", "판별", "개념",
            "맞는지", "옳은지", "틀린지",
        ],
        "patterns": [
            r"\b(?:true|false|t/f)\b",
            r"(?:맞는|옳은|틀린|참인|거짓)",
            r"다음\s*(?:중|문장).*(?:옳|맞|참|거짓|true|false)",
        ],
    },
    {
        "id": "frf",
        "label": "FRF 분해 (전달함수)",
        "menu": "2",
        "keywords": [
            "frf", "transfer function", "전달함수", "partial fraction",
            "부분분수", "ode", "bode", "impulse response",
            "g(s)", "h(s)", "frequency response",
            "omega_b", "omega_n", "break frequency",
            "subsystem", "서브시스템", "1차", "2차",
        ],
        "patterns": [
            r"x['\^]+.*=\s*f\(t\)",          # ODE form: x''' + ... = f(t)
            r"\d+\s*x['\^]",                  # coefficient + derivative
            r"s\^?\d+\s*[+\-]",              # polynomial in s
            r"전달\s*함수|transfer\s*function",
            r"부분\s*분수|partial\s*fraction",
            r"g\s*\(\s*s\s*\)|h\s*\(\s*s\s*\)",
        ],
    },
    {
        "id": "state_space",
        "label": "State-Space / 전이행렬 (상태방정식)",
        "menu": "3",
        "keywords": [
            "state space", "state-space", "상태공간", "상태방정식",
            "transition matrix", "전이행렬",
            "matrix exponential", "행렬지수", "expm",
            "state variable", "상태변수",
            "state equation", "first order", "1차",
            "phi(t)", "e^(at)", "convolution integral",
        ],
        "patterns": [
            r"state[\s-]?space|상태\s*공간|상태\s*방정식",
            r"transition\s*matrix|전이\s*행렬",
            r"matrix\s*exponential|행렬\s*지수",
            r"e\^?\s*\(\s*a\s*t\s*\)|expm|phi\s*\(\s*t\s*\)",
            r"ẋ\s*=\s*a\s*x|x\s*dot\s*=\s*a",
        ],
    },
    {
        "id": "inverse_laplace",
        "label": "Inverse Laplace Transform (역 라플라스)",
        "menu": "4",
        "keywords": [
            "inverse laplace", "laplace transform", "역 라플라스", "라플라스",
            "partial fraction", "부분분수",
            "f(s)", "transfer function", "전달함수",
            "poles", "극점", "residue", "유수",
            "heaviside", "unit step response",
        ],
        "patterns": [
            r"inverse\s*laplace|역\s*라플라스",
            r"laplace\s*transform|라플라스\s*변환",
            r"l\^?\s*[-−]?\s*1\s*\{",
            r"f\s*\(\s*s\s*\)\s*=",
            r"\d+\s*/\s*\(s",
        ],
    },
    {
        "id": "lagrange",
        "label": "Lagrange → 평형점 → 선형화 → 안정성",
        "menu": "5",
        "keywords": [
            "lagrange", "라그랑주", "라그랑지",
            "kinetic energy", "potential energy", "운동에너지", "위치에너지",
            "equilibrium", "평형점", "평형",
            "linearization", "선형화",
            "virtual work", "가상일",
            "generalized coordinate", "일반좌표", "일반화좌표",
        ],
        "patterns": [
            r"t\s*=.*(?:rdot|qdot|θ_dot|theta_dot)",  # T = ...qdot
            r"v\s*=.*(?:cos|sin|mg|kx)",               # V = potential
            r"lagrange|라그랑[주지]",
            r"평형\s*점|equilibrium\s*point",
            r"선형화|lineariz",
        ],
    },
    {
        "id": "modal",
        "label": "Modal Analysis (고유값/고유벡터)",
        "menu": "6",
        "keywords": [
            "modal", "모달", "eigenvalue", "고유값", "고유진동수",
            "eigenvector", "고유벡터", "mode shape", "모드형상",
            "natural frequency", "mass matrix", "stiffness matrix",
            "질량행렬", "강성행렬",
            "orthogonality", "직교성",
            "mass normalization", "질량정규화",
        ],
        "patterns": [
            r"k\s*u\s*=\s*[λλ]\s*m\s*u",
            r"modal\s*analysis|모달\s*분석|모달\s*해석",
            r"eigenvalue|고유값",
            r"mode\s*shape|모드\s*형상",
            r"m\s*=.*;\s*k\s*=",  # M = ...; K = ...
        ],
    },
    {
        "id": "damping",
        "label": "Proportional Damping Response",
        "menu": "7",
        "keywords": [
            "damping", "감쇠", "rayleigh", "레일리",
            "proportional", "비례감쇠",
            "alpha", "beta", "damping ratio", "감쇠비",
            "zeta", "nonproportional", "비비례",
            "damping matrix", "감쇠행렬",
        ],
        "patterns": [
            r"c\s*=\s*[αa].*m\s*\+\s*[βb].*k",  # C = αM + βK
            r"rayleigh\s*damping|레일리\s*감쇠|비례\s*감쇠",
            r"damping\s*ratio|감쇠\s*비",
            r"proportional\s*damp",
        ],
    },
    {
        "id": "base_excitation",
        "label": "Base Excitation (지반 가진)",
        "menu": "9",
        "keywords": [
            "base excitation", "지반가진", "지반 가진", "support motion",
            "ground motion", "transmissibility", "전달률",
            "isolation", "진동절연", "방진",
            "base motion", "moving support", "moving base",
            "y=asinwt", "y = a sin", "support excitation",
        ],
        "patterns": [
            r"base\s*excitat|지반\s*가진|support\s*(?:motion|excit)",
            r"transmissib|전달률",
            r"cy['\.]?\s*\+\s*ky|ky\s*\+\s*cy",
            r"y\s*=\s*[aA]\s*sin",
            r"ground\s*motion|moving\s*(?:base|support)",
            r"isolation|진동\s*절연|방진",
        ],
    },
    {
        "id": "convolution",
        "label": "Convolution Integral (Duhamel 적분)",
        "menu": "11",
        "keywords": [
            "convolution", "컨볼루��", "duhamel", "듀하멜",
            "impulse response", "충격응답", "임펄스응답",
            "piecewise", "구간별", "truncated",
            "underdamped", "미감쇠", "부족감쇠",
            "zero initial", "영초기조건",
        ],
        "patterns": [
            r"convolution\s*integral|duhamel",
            r"h\s*\(\s*t\s*[-−]\s*[τt]",
            r"∫.*h\s*\(.*\)\s*f\s*\(",
            r"impulse\s*response.*integr",
            r"piecewise|구간별",
            r"underdamp.*(?:response|system)",
        ],
    },
    {
        "id": "fourier",
        "label": "Fourier 급수 (주기 가진 응답)",
        "menu": "10",
        "keywords": [
            "fourier", "푸리에", "periodic", "주기", "square wave", "구형파",
            "sawtooth", "톱니파", "triangle wave", "삼각파",
            "periodic excitation", "주기가진", "주기 가진",
            "harmonic series", "조화급수",
            "first order", "first-order", "1차 시스템",
        ],
        "patterns": [
            r"fourier|푸리에",
            r"square[\s-]?wave|구형파",
            r"periodic\s*(?:excit|forc|load)|주기\s*(?:가진|하중)",
            r"c\s*x['\.]?\s*\+\s*k\s*x\s*=",
            r"f\s*\(\s*t\s*\)\s*=\s*f\s*\(\s*t\s*\+\s*t\s*\)",
            r"sawtooth|톱니|triangle\s*wave|삼각파",
        ],
    },
    {
        "id": "extended",
        "label": "확장 (Gyroscopic / Nonconservative)",
        "menu": "8",
        "keywords": [
            "gyroscopic", "자이로", "자이로스코픽",
            "nonconservative", "비보존", "circulatory", "순환력",
            "flutter", "divergence", "플러터", "발산",
            "skew-symmetric", "반대칭", "비대칭",
            "left eigenvector", "right eigenvector",
            "bi-orthogonality", "쌍직교",
            "kelvin-tait", "ziegler",
        ],
        "patterns": [
            r"gyroscop|자이로",
            r"nonconservat|비보존|circulat|순환",
            r"flutter|divergen|플러터|발산",
            r"skew[\s-]?symmetric|반대칭",
            r"bi[\s-]?orthogonal|쌍직교",
        ],
    },
]


def classify(text: str) -> RouteResult:
    """입력 텍스트를 분석하여 가장 적합한 문제 유형을 반환한다.

    Parameters
    ----------
    text : str
        문제 텍스트 또는 키워드.

    Returns
    -------
    RouteResult
        route, label, confidence, matched_keywords, suggestions.
    """
    text_lower = text.lower().strip()
    scores: list[tuple[str, str, float, list[str]]] = []

    for route in _ROUTES:
        score = 0.0
        matched: list[str] = []

        # Keyword matching
        for kw in route["keywords"]:
            if kw.lower() in text_lower:
                score += 2.0
                matched.append(kw)

        # Pattern matching (stronger signal)
        for pat in route.get("patterns", []):
            if re.search(pat, text_lower, re.IGNORECASE):
                score += 5.0
                matched.append(f"pattern:{pat[:30]}")

        scores.append((route["id"], route["label"], score, matched))

    scores.sort(key=lambda x: x[2], reverse=True)

    if not scores or scores[0][2] == 0:
        return RouteResult(
            route="unknown",
            label="분류 불가",
            confidence=0.0,
            suggestions=[s[0] for s in scores[:3] if s[2] > 0],
        )

    best_id, best_label, best_score, best_matched = scores[0]

    # Confidence: normalize by max possible score
    max_possible = max(best_score, 1.0)
    confidence = min(best_score / 20.0, 1.0)

    # Suggestions: other routes with nonzero score
    suggestions = [
        s[0] for s in scores[1:4] if s[2] > 0
    ]

    return RouteResult(
        route=best_id,
        label=best_label,
        confidence=confidence,
        matched_keywords=best_matched,
        suggestions=suggestions,
    )


def get_menu_number(route_id: str) -> str | None:
    """라우트 ID에 해당하는 메뉴 번호를 반환한다."""
    for route in _ROUTES:
        if route["id"] == route_id:
            return route.get("menu")
    return None


def list_routes() -> list[dict[str, str]]:
    """사용 가능한 라우트 목록을 반환한다."""
    return [{"id": r["id"], "label": r["label"], "menu": r["menu"]} for r in _ROUTES]
