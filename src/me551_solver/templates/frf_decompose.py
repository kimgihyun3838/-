"""FRF 분해 문제 템플릿 – ODE 문자열 파싱 + FRFSolver 위임.

Convenience wrapper that accepts either:
  1. An ODE string such as  "x''' + 3x'' + 6x' + 8x = f(t)"
  2. A raw coefficient list  [1, 3, 6, 8]

and delegates to ``core.frf.FRFSolver`` for the actual computation,
then packages the result with enhanced formatting.

ODE string grammar (accepted patterns)
---------------------------------------
The parser handles standard ME551 exam notation::

    a_n * x^(n) + ... + a_1 * x' + a_0 * x = f(t)

Supported derivative notations (all equivalent):
  * Prime marks:        x''', x'', x', x
  * d/dt notation:      d³x/dt³  (not required for exam problems)
  * Carets:             x^(3), x^(2), x^(1), x^(0)

Coefficient may be omitted (→ 1), or be a decimal/fraction.
RHS "= f(t)" (or "= F(t)" or "= u(t)") is ignored; only the LHS matters.

Examples
--------
>>> solver = FRFDecomposeSolver()
>>> result = solver.solve({"ode_string": "x''' + 3x'' + 6x' + 8x = f(t)"})
>>> result.final_answer["omega_1"]   # 2.0 rad/s
>>> result.final_answer["zeta_1"]    # 0.25

>>> result2 = solver.solve({"coefficients": [1, 1.2, 4.2, 4]})
>>> result2.final_answer["omega_1"]  # 2.0 rad/s
"""

from __future__ import annotations

import re
from typing import Any

from ..core.base import BaseSolver, SolverResult
from ..core.frf import FRFSolver

# ---------------------------------------------------------------------------
# ODE string parser
# ---------------------------------------------------------------------------

# Regex pieces
_COEFF   = r"([+-]?\s*(?:\d+(?:\.\d*)?|\.\d+)?)"   # optional signed number
_DERIV   = (
    r"x\s*(?:"
    r"'{1,9}"                       # prime marks: x', x'', x''', …
    r"|(?:\^\s*\(\s*(\d+)\s*\))"    # x^(n) notation
    r"|(?:\^\s*(\d+))"              # x^n  notation
    r")?"                           # no decoration → x^0
)


# Compiled pattern: optional coefficient, then the variable term
_TERM_RE = re.compile(
    r"([+-]?\s*(?:\d+(?:\.\d*)?)?)\s*"  # coefficient (possibly empty)
    r"\*?\s*"                            # optional explicit multiply
    r"(x\s*(?:'{1,9}|\^\s*\(\d+\)|\^\s*\d+)?)"  # x with decoration
    r"(?!\w)",                           # not followed by another word char
    re.IGNORECASE,
)


def _count_primes(token: str) -> int:
    """Return the derivative order encoded in a token like x'', x^(3), x^2."""
    # Prime marks
    primes = token.count("'")
    if primes:
        return primes
    # Caret notation:  x^(3) or x^3
    m = re.search(r"\^\s*\(?\s*(\d+)\s*\)?", token)
    if m:
        return int(m.group(1))
    # Plain x → 0th derivative
    return 0


def _parse_forcing_rhs(rhs: str) -> dict | None:
    """Parse the RHS of an ODE to extract forcing function parameters.

    Supported forms:
      - "f(t)", "F(t)", "u(t)"  → None (generic, no specific forcing)
      - "10sin(3t)", "10sin3t", "10sin(wt)", "10*sin(wt)"
      - "5cos(2t)", "5cos2t", "5*cos(2t)"
      - "Asin(wt)" where A is numeric, w is numeric or literal 'w'/'omega'

    Returns
    -------
    dict | None
        {"type": "harmonic", "F0": float, "omega": float|str, "func": "sin"|"cos"}
        or None if the RHS is generic / unparseable.
    """
    rhs = rhs.strip()
    if not rhs:
        return None

    # Skip generic forcing placeholders
    if re.match(r"^[fFuU]\s*\(\s*t\s*\)$", rhs):
        return None

    # Pattern: optional_coeff * sin/cos( omega * t )
    # Examples: 10sin(3t), 10*sin(3t), 10sin3t, sin(wt), 2.5cos(omega*t)
    m = re.match(
        r"([+-]?\s*(?:\d+(?:\.\d*)?|\.\d+)?)\s*\*?\s*"  # coefficient
        r"(sin|cos)\s*"                                    # trig function
        r"(?:\(\s*"                                        # optional opening paren
        r"([+-]?\d+(?:\.\d*)?|\.\d+|w|omega|ω)\s*\*?\s*t"  # omega * t
        r"\s*\)|"                                          # closing paren OR
        r"([+-]?\d+(?:\.\d*)?|\.\d+|w|omega|ω)\s*\*?\s*t)" # no parens
        r"\s*$",
        rhs,
        re.IGNORECASE,
    )
    if not m:
        return None

    raw_coeff = m.group(1).replace(" ", "")
    if raw_coeff in ("", "+"):
        F0 = 1.0
    elif raw_coeff == "-":
        F0 = -1.0
    else:
        F0 = float(raw_coeff)

    func = m.group(2).lower()
    raw_omega = (m.group(3) or m.group(4)).strip().lower()

    if raw_omega in ("w", "omega", "ω"):
        omega = raw_omega  # symbolic — caller decides
    else:
        omega = float(raw_omega)

    return {"type": "harmonic", "F0": F0, "omega": omega, "func": func}


def parse_ode_string(ode_string: str) -> tuple[list[float], dict | None]:
    """Parse an ODE string and return denominator coefficients and forcing info.

    The parser splits at '=', extracts coefficients from the LHS,
    and attempts to parse the RHS as a forcing function.

    Parameters
    ----------
    ode_string : str
        E.g. "x''' + 3x'' + 6x' + 8x = f(t)"
        E.g. "x'' + 0.6x' + 9x = 10sin(wt)"

    Returns
    -------
    tuple[list[float], dict | None]
        (coefficients [a_n, ..., a_0], forcing dict or None)

    Raises
    ------
    ValueError
        If no valid terms can be parsed.
    """
    # 1. Split LHS / RHS
    parts = ode_string.split("=")
    lhs = parts[0]
    rhs = parts[1].strip() if len(parts) > 1 else ""

    # 2. Normalise spacing around signs so tokenisation is easier
    lhs = re.sub(r"\s*([+-])\s*", r" \1", lhs).strip()
    # If the first character is a bare sign separate it for the regex
    lhs = re.sub(r"^\+", "", lhs)

    # 3. Find all terms
    matches = _TERM_RE.findall(lhs)
    if not matches:
        raise ValueError(
            f"Could not parse any x-terms from ODE string: '{ode_string}'"
        )

    coeff_map: dict[int, float] = {}
    for raw_coeff, raw_var in matches:
        raw_coeff = raw_coeff.replace(" ", "")
        # Parse coefficient
        if raw_coeff in ("", "+"):
            c = 1.0
        elif raw_coeff == "-":
            c = -1.0
        else:
            c = float(raw_coeff)

        order = _count_primes(raw_var.replace(" ", ""))
        coeff_map[order] = coeff_map.get(order, 0.0) + c

    if not coeff_map:
        raise ValueError(
            f"No terms with variable 'x' found in ODE string: '{ode_string}'"
        )

    # 4. Build coefficient list a_n, a_{n-1}, ..., a_0
    max_order = max(coeff_map.keys())
    coefficients = [coeff_map.get(n, 0.0) for n in range(max_order, -1, -1)]

    # 5. Parse RHS forcing function
    forcing = _parse_forcing_rhs(rhs) if rhs else None

    return coefficients, forcing


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------


class FRFDecomposeSolver(BaseSolver):
    """FRF 분해 편의 래퍼 – ODE 문자열 또는 계수 리스트를 받아 FRFSolver에 위임.

    Provides:
      * ``ode_string`` → parsed coefficient list
      * ``coefficients`` → forwarded directly
      * Enhanced result with flattened first-order / second-order parameters
    """

    def solve(self, params: dict) -> SolverResult:
        """Perform FRF decomposition.

        Parameters
        ----------
        params : dict
            ode_string   : str | None
                ODE in human-readable notation, e.g.
                "x''' + 3x'' + 6x' + 8x = f(t)"
            coefficients : list[float] | None
                Denominator polynomial coefficients [a_n, …, a_0].
                Required if ``ode_string`` is not provided.
            numerator    : list[float] | None
                Numerator coefficients (default [1]).
            omega_range  : tuple | None
                (omega_min, omega_max[, num_points]) for numeric G(jω).
        """
        ode_str     = params.get("ode_string")
        coefficients = params.get("coefficients")
        numerator   = params.get("numerator", [1])
        omega_range = params.get("omega_range")

        steps: list[tuple[str, Any]] = []
        given: dict[str, Any] = {}

        # ---------------------------------------------------------------
        # 1. Resolve coefficients
        # ---------------------------------------------------------------
        parsed_forcing = None
        if ode_str is not None:
            given["ode_string"] = ode_str
            parsed_coeffs, parsed_forcing = parse_ode_string(ode_str)
            given["coefficients"] = parsed_coeffs
            parse_detail = (
                f"Input:   '{ode_str}'\n"
                f"Coefficients [a_n, …, a_0]: {parsed_coeffs}"
            )
            if parsed_forcing:
                parse_detail += (
                    f"\nForcing function detected: "
                    f"{parsed_forcing['F0']} {parsed_forcing['func']}"
                    f"({parsed_forcing['omega']}·t)"
                )
            steps.append(("ODE string parsed", parse_detail))
        elif coefficients is not None:
            parsed_coeffs = [float(c) for c in coefficients]
            given["coefficients"] = parsed_coeffs
            steps.append((
                "Coefficients provided directly",
                f"[a_n, …, a_0] = {parsed_coeffs}"
            ))
        else:
            raise ValueError(
                "Provide either 'ode_string' or 'coefficients'."
            )

        n = len(parsed_coeffs) - 1
        given["system_order"] = n
        given["numerator"] = numerator

        # ---------------------------------------------------------------
        # 2. Delegate to FRFSolver
        # ---------------------------------------------------------------
        frf_solver = FRFSolver()
        frf_params: dict[str, Any] = {
            "coefficients": parsed_coeffs,
            "numerator":    numerator,
        }
        if omega_range is not None:
            frf_params["omega_range"] = omega_range

        # Use explicitly provided forcing, or auto-parsed from ODE string
        forcing = params.get("forcing")
        if forcing is None and parsed_forcing is not None:
            forcing = parsed_forcing
        if forcing is not None:
            frf_params["forcing"] = forcing

        base_result = frf_solver.solve(frf_params)

        # Merge FRF steps into our steps (with clear labels)
        for desc, expr in base_result.steps:
            steps.append((desc, expr))

        # ---------------------------------------------------------------
        # 3. Flatten and enhance final_answer
        # ---------------------------------------------------------------
        fa = dict(base_result.final_answer)   # copy

        # Convenience aliases for exam-style answers
        first_order  = fa.get("first_order",  [])
        second_order = fa.get("second_order", [])

        if first_order:
            fo = first_order[0]
            fa["pole_real"]   = -fo["omega_b"]
            fa["bandwidth_G1"] = fo["omega_b"]

        if second_order:
            so = second_order[0]
            fa["omega_n"]  = so["omega_n"]
            fa["zeta"]     = so["zeta"]
            fa["omega_1"]  = so["omega_n"]  # alias
            fa["zeta_1"]   = so["zeta"]     # alias

        # DC gain shortcut
        fa["G_dc"] = fa.get("dc_gain", None)

        # ---------------------------------------------------------------
        # 4. Summary step
        # ---------------------------------------------------------------
        summary_lines = [
            f"System order: {n}",
            f"Characteristic polynomial: {fa.get('transfer_function', '')}",
        ]
        if first_order:
            fo = first_order[0]
            summary_lines.append(
                f"G₁(s): 1st-order,  ω_b = {fo['omega_b']:.6g} rad/s"
            )
        if second_order:
            so = second_order[0]
            summary_lines.append(
                f"G₂(s): 2nd-order,  ω_n = {so['omega_n']:.6g} rad/s,  "
                f"ζ = {so['zeta']:.6g}"
            )
        steps.append(("Summary", "\n".join(summary_lines)))

        # ---------------------------------------------------------------
        # 5. Build and return SolverResult
        # ---------------------------------------------------------------
        return SolverResult(
            problem_type="FRF Decomposition (ODE → Transfer Function → Partial Fractions)",
            given=given,
            steps=steps,
            final_answer=fa,
            sanity_check=base_result.sanity_check,
        )

    def get_input_template(self) -> dict:
        """Return the expected input parameter schema."""
        return {
            "ode_string": {
                "type": "str | None",
                "description": (
                    "ODE in human-readable notation.\n"
                    "Examples:\n"
                    "  'x''' + 3x'' + 6x' + 8x = f(t)'\n"
                    "  'x'' + 0.2x' + 4x = f(t)'\n"
                    "Either this or 'coefficients' must be supplied."
                ),
                "example": "x''' + 3x'' + 6x' + 8x = f(t)",
            },
            "coefficients": {
                "type": "list[float] | None",
                "description": (
                    "Denominator polynomial coefficients [a_n, …, a_0], "
                    "high-to-low order.\n"
                    "Example for x'''+3x''+6x'+8x: [1, 3, 6, 8].\n"
                    "Either this or 'ode_string' must be supplied."
                ),
                "example": [1, 3, 6, 8],
            },
            "numerator": {
                "type": "list[float] | None",
                "description": "Numerator polynomial coefficients (default [1]).",
                "default": [1],
            },
            "omega_range": {
                "type": "tuple | None",
                "description": (
                    "(omega_min, omega_max[, num_points]) to compute |G(jω)| "
                    "and phase numerically."
                ),
                "example": [0.01, 100, 500],
            },
        }
