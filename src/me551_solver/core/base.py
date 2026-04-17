"""Base interfaces for all ME551 solvers."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any


@dataclass
class SolverResult:
    """모든 솔버가 반환하는 통일된 결과 구조체.

    Attributes:
        problem_type: 문제 유형 (예: "FRF", "Lagrange", "Modal Analysis")
        given: 주어진 조건 dict
        steps: 풀이 과정 리스트 – 각 원소는 (설명, 수식/값) 튜플
        final_answer: 최종 답 dict
        sanity_check: 차원 검증·극한 검증 등 결과 문자열
    """

    problem_type: str
    given: dict[str, Any] = field(default_factory=dict)
    steps: list[tuple[str, Any]] = field(default_factory=list)
    final_answer: dict[str, Any] = field(default_factory=dict)
    sanity_check: str = ""


class BaseSolver(ABC):
    """모든 문제 유형 솔버의 추상 베이스 클래스."""

    @abstractmethod
    def solve(self, params: dict) -> SolverResult:
        """주어진 파라미터로 문제를 풀고 SolverResult를 반환한다."""
        pass

    @abstractmethod
    def get_input_template(self) -> dict:
        """필요한 입력 파라미터 설명 dict를 반환한다."""
        pass
