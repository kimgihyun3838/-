# ConceptDB — T/F Verdict

## Given
- **statement**: Linear systems hold principle of superposition, i.e., both homogeneity andadditivity, but nonlinear systems do not hold at least one of homogeneity and additivity

## Derivation Steps
### Step 1: Query tokens

```
['additivity', 'and', 'andadditivity', 'both', 'but', 'hold', 'homogeneity', 'least', 'linear', 'nonlinear', 'not', 'one', 'principle', 'superposition', 'systems']
```

### Step 2: Match #1: rule_001 (score=33.0)

```
{'statement': 'A linear system satisfies both homogeneity and additivity (superposition principle).', 'topic': 'linear_systems', 'keyword_overlap': ['additivity', 'homogeneity', 'superposition']}
```

### Step 3: Match #2: rule_005 (score=23.0)

```
{'statement': 'The principle of superposition applies to nonlinear systems when the excitation amplitude is sufficiently small.', 'topic': 'linear_systems', 'keyword_overlap': ['nonlinear', 'superposition']}
```

### Step 4: Match #3: rule_002 (score=22.0)

```
{'statement': 'If a system satisfies homogeneity alone, it is necessarily linear.', 'topic': 'linear_systems', 'keyword_overlap': ['additivity', 'homogeneity', 'superposition']}
```

## Final Answer
- **rule_id**: `rule_001`
- **topic**: `linear_systems`
- **matched_statement**: `A linear system satisfies both homogeneity and additivity (superposition principle).`
- **verdict**: `TRUE`
- **reason**: `Linearity is defined by two properties: homogeneity (scaling input scales output by the same factor) and additivity (sum of inputs produces sum of individual outputs). Together these constitute the superposition principle. All systems described by linear ODEs with constant coefficients satisfy both. A single property alone is insufficient to guarantee linearity.`
- **lecture_ref**: `Lecture 2`
- **match_score**: `33.0`
- **common_trap**: `Students sometimes confuse 'superposition' as a single condition; it is actually the combination of both homogeneity and additivity.`
- **related_rules**: `['rule_002', 'rule_003']`
- **also_consider**: `[{'rule_id': 'rule_005', 'statement': 'The principle of superposition applies to nonlinear systems when the excitation amplitude is sufficiently small.', 'verdict': 'FALSE', 'topic': 'linear_systems'}, {'rule_id': 'rule_002', 'statement': 'If a system satisfies homogeneity alone, it is necessarily linear.', 'verdict': 'FALSE', 'topic': 'linear_systems'}]`

## Sanity Check
Confidence: HIGH — strong keyword overlap with the matched rule. (score=33.0, 60 rules loaded)
