"""T/F 개념 판별 솔버 v2 — 구조적 문장 분석 + 충돌 검사 기반.

v1: 키워드 매칭만으로 가장 유사한 rule의 verdict 반환
v2: 문장 구조 파싱(부정어/조건/수치/범위) + 규칙과의 충돌 검사
    → 함정 문장 감지 가능, Phase 1 advisory 모드 지원

Loads rules from data/tf_rules.yaml.  Fully offline — requires only pyyaml.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml

from .base import BaseSolver, SolverResult
from .paths import get_data_file

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
_TF_RULES_PATH = Path(get_data_file("tf_rules.yaml"))

# ---------------------------------------------------------------------------
# Statement parsing: negation, scope, conditions, numerics
# ---------------------------------------------------------------------------

_NEGATION_WORDS = frozenset({
    "not", "cannot", "never", "no", "doesn't", "don't",
    "isn't", "aren't", "won't", "unable", "impossible",
    "neither", "nor", "can't", "shouldn't", "couldn't",
})

_NEGATION_PHRASES = [
    "does not", "do not", "can not", "cannot", "is not", "are not",
    "will not", "should not", "need not", "cannot be", "is never",
    "not applicable", "not valid", "not possible", "does no",
]

_SCOPE_UNIVERSAL = frozenset({
    "always", "all", "every", "any", "both", "entire", "whole",
    "general", "arbitrary", "regardless",
})

_SCOPE_RESTRICTIVE = frozenset({
    "only", "solely", "exclusively", "just", "merely",
})


def _parse_statement(text: str) -> dict[str, Any]:
    """Parse a T/F statement into structural components for comparison."""
    text_lower = text.lower().strip()

    # --- 0. Strip contrastive clauses before negation detection ---
    # Sentences like "X holds, but Y does not" use negation only to
    # contrast; the overall statement is affirmative.  Remove the
    # contrastive tail so that "do not" in "but ... do not" is not
    # counted as a negation of the main claim.
    main_clause = text_lower
    contrastive_patterns = [
        r",?\s+but\s+.+$",          # ", but nonlinear systems do not"
        r",?\s+whereas\s+.+$",      # ", whereas nonlinear ..."
        r",?\s+while\s+.+$",        # ", while nonlinear ..."
        r",?\s+however\s+.+$",      # ", however nonlinear ..."
    ]
    for cp in contrastive_patterns:
        main_clause = re.sub(cp, "", main_clause)

    # --- 1. Negation (on main clause only) ---
    found_negations: list[str] = []
    for phrase in _NEGATION_PHRASES:
        if phrase in main_clause:
            found_negations.append(phrase)
    # Individual negation words (avoid double-counting from phrases)
    phrase_words = set()
    for p in found_negations:
        phrase_words.update(p.split())
    for word in _NEGATION_WORDS:
        if re.search(r"\b" + re.escape(word) + r"\b", main_clause):
            if word not in phrase_words:
                found_negations.append(word)

    # --- 2. Scope ---
    scope_universal: list[str] = []
    scope_restrictive: list[str] = []
    for word in _SCOPE_UNIVERSAL:
        if re.search(r"\b" + re.escape(word) + r"\b", text_lower):
            scope_universal.append(word)
    for word in _SCOPE_RESTRICTIVE:
        if re.search(r"\b" + re.escape(word) + r"\b", text_lower):
            scope_restrictive.append(word)

    # --- 3. Conditions ---
    conditions: list[str] = []
    cond_patterns = [
        r"\bfor\s+(free|forced|harmonic|transient|arbitrary|steady[\s-]?state"
        r"|undamped|damped|nonlinear|linear|all|any|both)[^,;.]*",
        r"\bunder\s+(harmonic|arbitrary|periodic|transient|forced|free)[^,;.]*",
        r"\bduring\s+(free|forced|transient|steady)[^,;.]*",
        r"\bin\s+the\s+(time|frequency)\s+domain",
        r"\bin\s+(free|forced|transient|steady)\s+[^,;.]*",
        r"\bincluding\s+(free|forced|transient)[^,;.]*",
    ]
    for pat in cond_patterns:
        for m in re.finditer(pat, text_lower):
            conditions.append(m.group(0).strip())

    # --- 4. Dimension / numeric patterns ---
    dimensions: dict[str, str] = {}
    dim_patterns = [
        (r"\bp\s*[*×·]\s*n\b", "p*n"),
        (r"\bpn\b", "p*n"),
        (r"\b2\s*[*×·]?\s*n\b", "2n"),
        (r"\b2n\b", "2n"),
        (r"\bn\s*[*×·]\s*p\b", "p*n"),
        (r"\bp\b(?=[\s-]*(?:dimensional|dimension|차원))", "p"),
        (r"\bn\b(?=[\s-]*(?:dimensional|dimension|차원))", "n"),
        (r"\b2n\b(?=[\s-]*(?:dimensional|dimension|차원))", "2n"),
        (r"\bp\s*[*×·]?\s*n\b(?=[\s-]*(?:dimensional|dimension|차원))", "p*n"),
    ]
    for pat, label in dim_patterns:
        if re.search(pat, text_lower):
            dimensions["state_dimension"] = label

    # --- 5. Domain keywords ---
    domain_keywords: set[str] = set()
    domain_map = {
        "time domain": "time_domain",
        "frequency domain": "frequency_domain",
        "convolution": "convolution",
        "multiplication": "multiplication",
        "free vibration": "free_vibration",
        "free response": "free_vibration",
        "forced vibration": "forced_vibration",
        "forced response": "forced_vibration",
        "transient": "transient",
        "steady state": "steady_state",
        "steady-state": "steady_state",
        "harmonic": "harmonic",
    }
    for term, key in domain_map.items():
        if term in text_lower:
            domain_keywords.add(key)

    return {
        "has_negation": len(found_negations) > 0,
        "negation_count": len(found_negations),
        "negation_words": found_negations,
        "scope_universal": scope_universal,
        "scope_restrictive": scope_restrictive,
        "conditions": conditions,
        "dimensions": dimensions,
        "domain_keywords": domain_keywords,
        "text_lower": text_lower,
    }


def _detect_conflicts(
    input_parsed: dict[str, Any],
    rule: dict[str, Any],
    rule_parsed: dict[str, Any],
) -> list[dict[str, Any]]:
    """Detect structural conflicts between input statement and matched rule.

    Returns a list of conflict dicts, each with type, detail, and severity.
    """
    conflicts: list[dict[str, Any]] = []
    input_text = input_parsed["text_lower"]

    # ----- 1. Negation mismatch -----
    input_neg = input_parsed["has_negation"]
    rule_neg = rule_parsed["has_negation"]
    if input_neg != rule_neg:
        conflicts.append({
            "type": "negation",
            "detail": (
                f"입력 문장: 부정어 {'있음' if input_neg else '없음'} "
                f"{input_parsed['negation_words']}, "
                f"규칙 문장: 부정어 {'있음' if rule_neg else '없음'} "
                f"{rule_parsed['negation_words']}"
            ),
            "severity": "HIGH",
        })

    # ----- 2. Trap patterns (rule-defined regex traps) -----
    trap_patterns = rule.get("trap_patterns", [])
    for trap in trap_patterns:
        pattern = trap.get("pattern", "")
        if pattern and re.search(pattern, input_text):
            conflicts.append({
                "type": "trap_pattern",
                "detail": trap.get("explanation", f"함정 패턴 매칭: {pattern}"),
                "severity": "HIGH",
            })

    # ----- 3. Condition conflict (invalid_for) -----
    invalid_conditions = rule.get("invalid_for", [])
    for inv_cond in invalid_conditions:
        inv_lower = inv_cond.lower()
        if inv_lower in input_text:
            conflicts.append({
                "type": "condition",
                "detail": f"입력에 '{inv_cond}' 포함 — 규칙의 적용 제외 조건",
                "severity": "HIGH",
            })

    # ----- 4. Valid_for scope check -----
    valid_conditions = rule.get("valid_for", [])
    if valid_conditions:
        # Check for domain keywords that contradict the valid scope
        contradictions = [
            ("free_vibration", ["forced vibration", "forced response", "forced"]),
            ("transient", ["harmonic", "steady state", "steady-state"]),
            ("convolution", ["multiplication"]),
        ]
        for domain_kw in input_parsed.get("domain_keywords", set()):
            for bad_kw, valid_scopes in contradictions:
                if domain_kw == bad_kw:
                    for vs in valid_scopes:
                        if any(vs.lower() in vc.lower() for vc in valid_conditions):
                            conflicts.append({
                                "type": "scope_conflict",
                                "detail": (
                                    f"입력: '{domain_kw.replace('_', ' ')}' 언급, "
                                    f"규칙 유효 범위: {', '.join(valid_conditions)}"
                                ),
                                "severity": "HIGH",
                            })
                            break

    # ----- 5. Numerical/dimension mismatch -----
    key_quantities = rule.get("key_quantities", {})
    input_dims = input_parsed.get("dimensions", {})
    if key_quantities and input_dims:
        for qty_name, correct_val in key_quantities.items():
            input_val = input_dims.get(qty_name)
            if input_val and input_val != correct_val:
                conflicts.append({
                    "type": "numerical",
                    "detail": f"규칙: {qty_name}={correct_val}, 입력: {qty_name}={input_val}",
                    "severity": "HIGH",
                })

    # ----- 6. Scope overextension -----
    if valid_conditions and input_parsed.get("scope_universal"):
        overextend_words = {"both", "always", "any", "arbitrary", "regardless", "all"}
        used = set(input_parsed["scope_universal"]) & overextend_words
        if used:
            conflicts.append({
                "type": "scope_overextension",
                "detail": (
                    f"범위 확장 표현: {', '.join(used)} — "
                    f"규칙 유효 범위: {', '.join(valid_conditions)}"
                ),
                "severity": "MEDIUM",
            })

    return conflicts


# ---------------------------------------------------------------------------
# Keyword tokenizer (unchanged from v1)
# ---------------------------------------------------------------------------

_STOPWORDS = frozenset({
    # English function words that pollute matching
    "the", "and", "for", "are", "was", "has", "its", "can", "not", "but",
    "with", "from", "this", "that", "than", "then", "they", "them", "their",
    "there", "these", "those", "been", "being", "have", "had", "does", "did",
    "will", "would", "could", "should", "shall", "may", "might", "must",
    "also", "each", "every", "both", "such", "when", "where", "which",
    "while", "about", "into", "through", "during", "before", "after",
    "above", "below", "between", "under", "over", "same", "other",
    "more", "most", "some", "any", "all", "only", "own", "very",
    "just", "because", "here", "how", "what", "who", "whom",
    "give", "gives", "given", "used", "using", "use",
    "says", "said", "called", "known", "means", "following",
})


def _tokenize(text: str) -> set[str]:
    """Lower-case words (>=3 chars, no stopwords) from *text* for matching."""
    if not text:
        return set()
    tokens = re.findall(r"[a-zA-Z가-힣]{3,}", text)
    return {t.lower() for t in tokens} - _STOPWORDS


# ---------------------------------------------------------------------------
# ConceptDBSolver
# ---------------------------------------------------------------------------

class ConceptDBSolver(BaseSolver):
    """강의 노트 기반 T/F 개념 판별 솔버 v2.

    v2 주요 변경:
    - 문장 구조 분석 (부정어, 조건, 수치, 범위)
    - 충돌 검사 (입력 vs 매칭 규칙 비교)
    - Phase 1 advisory: 개념 + 주의 포인트 + 추정 판정
    """

    def __init__(self) -> None:
        self._rules: list[dict[str, Any]] = []
        self._rules_by_id: dict[str, dict] = {}
        self._topics: set[str] = set()
        self._load_rules()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _load_rules(self) -> None:
        """Load tf_rules.yaml into memory."""
        if not _TF_RULES_PATH.exists():
            raise FileNotFoundError(
                f"T/F rules database not found at: {_TF_RULES_PATH}\n"
                "Ensure data/tf_rules.yaml is present in the project."
            )
        with _TF_RULES_PATH.open(encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
        rules_list: list[dict] = raw.get("rules", [])
        if not rules_list:
            raise ValueError("tf_rules.yaml loaded but contained no rules.")
        for rule in rules_list:
            rule_id = rule.get("id", "")
            self._rules.append(rule)
            if rule_id:
                self._rules_by_id[rule_id] = rule
            topic = rule.get("topic", "")
            if topic:
                self._topics.add(topic)

    # ------------------------------------------------------------------
    # BaseSolver interface
    # ------------------------------------------------------------------

    def solve(self, params: dict) -> SolverResult:
        """Match a statement to the best-fitting T/F rule with conflict analysis.

        Parameters
        ----------
        params : dict
            - ``"statement"`` (str, required): concept statement to evaluate
            - ``"keywords"`` (list[str], optional): extra boost keywords
            - ``"topic"`` (str, optional): restrict to specific topic

        Returns
        -------
        SolverResult with conflict analysis and advisory output
        """
        statement: str = params.get("statement", "")
        extra_keywords: list[str] = params.get("keywords", [])
        topic_filter: str | None = params.get("topic", None)

        if not statement:
            return SolverResult(
                problem_type="ConceptDB",
                given=params,
                steps=[("Error", "No statement provided.")],
                final_answer={"error": "statement is required"},
                sanity_check="No query given.",
            )

        # Parse input statement structure
        input_parsed = _parse_statement(statement)

        # Build search tokens
        query_tokens = _tokenize(statement)
        for kw in extra_keywords:
            query_tokens.update(_tokenize(kw))

        # Candidate pool
        candidates = self._rules
        if topic_filter:
            candidates = [r for r in candidates if r.get("topic", "") == topic_filter]
            if not candidates:
                candidates = [
                    r for r in self._rules
                    if topic_filter.lower() in r.get("topic", "").lower()
                ]
        if not candidates:
            return SolverResult(
                problem_type="ConceptDB",
                given=params,
                steps=[("Info", f"No rules found for topic '{topic_filter}'.")],
                final_answer={"error": f"No rules in topic '{topic_filter}'"},
                sanity_check="Topic filter yielded no candidates.",
            )

        # Score each rule
        scored = self._score_rules(query_tokens, candidates)

        # Build result with structural conflict analysis
        return self._build_result(params, scored, query_tokens, input_parsed)

    def get_input_template(self) -> dict:
        return {
            "statement": "판별할 개념 문장 (str) — 영어 또는 한국어",
            "keywords": "(선택) 추가 키워드 리스트 (list[str])",
            "topic": (
                "(선택) 토픽 필터 (str). 가능한 값: "
                + ", ".join(sorted(self._topics))
            ),
        }

    # ------------------------------------------------------------------
    # Public search helpers
    # ------------------------------------------------------------------

    def search_by_keywords(self, keywords: list[str]) -> list[dict]:
        """Return rules ranked by keyword overlap with *keywords*."""
        query_tokens: set[str] = set()
        for kw in keywords:
            query_tokens.update(_tokenize(kw))
        scored = self._score_rules(query_tokens, self._rules)
        return [r for r, _ in scored if _ > 0]

    def search_by_topic(self, topic: str) -> list[dict]:
        """Return all rules for *topic* (exact or partial match)."""
        exact = [r for r in self._rules if r.get("topic", "") == topic]
        if exact:
            return exact
        return [r for r in self._rules if topic.lower() in r.get("topic", "").lower()]

    def get_rule_by_id(self, rule_id: str) -> dict | None:
        """Return a specific rule by its id (e.g. ``'rule_036'``)."""
        return self._rules_by_id.get(rule_id)

    def list_topics(self) -> list[str]:
        """Return sorted list of all available topics."""
        return sorted(self._topics)

    def list_all_rules(self) -> list[dict]:
        """Return all loaded rules (read-only snapshot)."""
        return list(self._rules)

    # ------------------------------------------------------------------
    # Internal: scoring (same weights as v1)
    # ------------------------------------------------------------------

    def _score_rules(
        self,
        query_tokens: set[str],
        candidates: list[dict],
    ) -> list[tuple[dict, float]]:
        """Score each candidate rule by keyword overlap.

        Weights: keyword ×3, statement ×2, topic ×1, reason ×1.
        """
        results: list[tuple[dict, float]] = []
        for rule in candidates:
            score = 0.0
            rule_keywords = {kw.lower() for kw in rule.get("keywords", [])}
            score += len(query_tokens & rule_keywords) * 3.0
            stmt_tokens = _tokenize(rule.get("statement", ""))
            score += len(query_tokens & stmt_tokens) * 2.0
            topic_tokens = _tokenize(rule.get("topic", ""))
            score += len(query_tokens & topic_tokens) * 1.0
            reason_tokens = _tokenize(rule.get("reason", ""))
            score += len(query_tokens & reason_tokens) * 1.0
            results.append((rule, score))
        results.sort(key=lambda x: x[1], reverse=True)
        return results

    # ------------------------------------------------------------------
    # Internal: result builder with conflict analysis (v2)
    # ------------------------------------------------------------------

    def _build_result(
        self,
        params: dict,
        scored: list[tuple[dict, float]],
        query_tokens: set[str],
        input_parsed: dict[str, Any],
    ) -> SolverResult:
        """Convert scored rule list into SolverResult with conflict analysis."""

        if not scored or scored[0][1] == 0.0:
            return SolverResult(
                problem_type="ConceptDB -- T/F 분석",
                given=params,
                steps=[("Warning", "No rules matched the given statement/keywords.")],
                final_answer={
                    "verdict": None,
                    "message": "매칭 규칙 없음. 키워드를 변경하거나 토픽으로 검색하세요.",
                    "available_topics": self.list_topics(),
                },
                sanity_check="No matches.",
            )

        best_rule, best_score = scored[0]

        # Parse the matched rule's statement for structural comparison
        rule_parsed = _parse_statement(best_rule.get("statement", ""))

        # Detect conflicts between input and matched rule
        conflicts = _detect_conflicts(input_parsed, best_rule, rule_parsed)

        # Determine verdict
        rule_verdict = best_rule.get("verdict")
        rule_verdict_str = (
            "TRUE" if rule_verdict is True
            else "FALSE" if rule_verdict is False
            else "UNKNOWN"
        )

        # Apply conflict-based adjustment
        high_conflicts = [c for c in conflicts if c["severity"] == "HIGH"]
        has_high_conflict = len(high_conflicts) > 0
        has_any_conflict = len(conflicts) > 0

        if has_high_conflict:
            # If the matched rule is already FALSE and the conflicts are
            # only trap_pattern matches, the traps *confirm* that the input
            # is indeed the false variant → keep FALSE (do not flip).
            trap_only = all(c["type"] == "trap_pattern" for c in high_conflicts)
            if rule_verdict is False and trap_only:
                adjusted_verdict = "FALSE"
                confidence = "CONFLICT_DETECTED"
                verdict_note = (
                    f"매칭 규칙({best_rule['id']})은 {rule_verdict_str}이며, "
                    f"함정 패턴이 {len(high_conflicts)}개 매칭됨 "
                    f"-> 판정 유지: {adjusted_verdict}"
                )
            else:
                adjusted_verdict = "FALSE" if rule_verdict is True else "TRUE"
                confidence = "CONFLICT_DETECTED"
                verdict_note = (
                    f"매칭 규칙({best_rule['id']})은 {rule_verdict_str}이지만, "
                    f"입력 문장과 {len(high_conflicts)}개의 구조적 차이가 감지됨 "
                    f"-> 추정: {adjusted_verdict}"
                )
        elif has_any_conflict:
            adjusted_verdict = rule_verdict_str
            confidence = "CAUTION"
            verdict_note = (
                f"규칙({best_rule['id']}) verdict: {rule_verdict_str}, "
                f"일부 차이 감지 — 주의 필요"
            )
        else:
            adjusted_verdict = rule_verdict_str
            if best_score >= 9:
                confidence = "HIGH"
            elif best_score >= 4:
                confidence = "MEDIUM"
            else:
                confidence = "LOW"
            verdict_note = f"규칙({best_rule['id']})과 구조 일치"

        # ---- Build steps ----
        steps: list[tuple[str, Any]] = []
        steps.append(("Query tokens", sorted(query_tokens)))

        # Input structure analysis
        struct_lines: list[str] = []
        if input_parsed["has_negation"]:
            struct_lines.append(f"부정어: {input_parsed['negation_words']}")
        if input_parsed["scope_universal"]:
            struct_lines.append(f"범위 확장: {input_parsed['scope_universal']}")
        if input_parsed["scope_restrictive"]:
            struct_lines.append(f"범위 제한: {input_parsed['scope_restrictive']}")
        if input_parsed["conditions"]:
            struct_lines.append(f"조건절: {input_parsed['conditions']}")
        if input_parsed["dimensions"]:
            struct_lines.append(f"수치/차원: {input_parsed['dimensions']}")
        if input_parsed["domain_keywords"]:
            struct_lines.append(
                f"도메인: {', '.join(sorted(input_parsed['domain_keywords']))}"
            )
        if struct_lines:
            steps.append(("입력 문장 구조 분석", "\n".join(struct_lines)))

        # Top matches
        for idx, (rule, score) in enumerate(scored[:3]):
            if score <= 0:
                break
            rule_kws = {kw.lower() for kw in rule.get("keywords", [])}
            overlap = sorted(query_tokens & rule_kws)
            steps.append((
                f"Match #{idx + 1}: {rule['id']} (score={score:.1f})",
                {
                    "statement": rule.get("statement", ""),
                    "topic": rule.get("topic", ""),
                    "keyword_overlap": overlap,
                },
            ))

        # Conflict analysis
        if conflicts:
            conflict_lines = []
            for c in conflicts:
                conflict_lines.append(f"[{c['severity']}] {c['type']}: {c['detail']}")
            steps.append(("충돌 분석", "\n".join(conflict_lines)))

        # ---- Final answer ----
        final_answer: dict[str, Any] = {
            "rule_id": best_rule.get("id"),
            "topic": best_rule.get("topic"),
            "matched_statement": best_rule.get("statement"),
            "rule_verdict": rule_verdict_str,
            "adjusted_verdict": adjusted_verdict,
            "confidence": confidence,
            "verdict_note": verdict_note,
            "reason": best_rule.get("reason", "").strip(),
            "lecture_ref": best_rule.get("lecture_ref", ""),
            "match_score": best_score,
            "conflicts": conflicts,
        }

        if best_rule.get("common_trap"):
            final_answer["common_trap"] = best_rule["common_trap"]
        if best_rule.get("related_rules"):
            final_answer["related_rules"] = best_rule["related_rules"]
        if best_rule.get("valid_for"):
            final_answer["valid_for"] = best_rule["valid_for"]
        if best_rule.get("invalid_for"):
            final_answer["invalid_for"] = best_rule["invalid_for"]

        # Alternate matches
        top_n = [r for r, s in scored[:3] if s > 0]
        if len(top_n) > 1:
            final_answer["also_consider"] = [
                {
                    "rule_id": r.get("id"),
                    "statement": r.get("statement", ""),
                    "verdict": "TRUE" if r.get("verdict") is True else "FALSE",
                    "topic": r.get("topic", ""),
                }
                for r in top_n[1:]
            ]

        # Advisory block (Phase 1 output)
        advisory_lines: list[str] = []
        advisory_lines.append(f"매칭 규칙: {best_rule.get('statement', '')}")
        advisory_lines.append(f"규칙 verdict: {rule_verdict_str}")
        if conflicts:
            advisory_lines.append("주의 포인트:")
            for c in conflicts:
                advisory_lines.append(f"  - [{c['type']}] {c['detail']}")
            if has_high_conflict:
                advisory_lines.append(f"추정 판정: {adjusted_verdict} (규칙과 반대)")
        else:
            advisory_lines.append(f"추정 판정: {adjusted_verdict}")
        if best_rule.get("common_trap"):
            advisory_lines.append(f"흔한 함정: {best_rule['common_trap']}")
        if best_rule.get("valid_for"):
            advisory_lines.append(f"유효 조건: {', '.join(best_rule['valid_for'])}")
        if best_rule.get("invalid_for"):
            advisory_lines.append(f"무효 조건: {', '.join(best_rule['invalid_for'])}")

        final_answer["advisory"] = "\n".join(advisory_lines)

        # Backward compatibility: keep "verdict" field (now uses adjusted)
        final_answer["verdict"] = adjusted_verdict

        sanity = (
            f"Confidence: {confidence} | Score: {best_score:.1f} | "
            f"Conflicts: {len(conflicts)} ({len(high_conflicts)} HIGH) | "
            f"Rules: {len(self._rules)}"
        )

        return SolverResult(
            problem_type="ConceptDB -- T/F 분석",
            given=params,
            steps=steps,
            final_answer=final_answer,
            sanity_check=sanity,
        )
