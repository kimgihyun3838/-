"""안정성 분석 솔버 -- 고유값 기반 안정성 판별.

Supports two modes:
  1. "linearized_ode" -- from linearized EOM coefficients (mass, damping, stiffness)
  2. "state_matrix"   -- from a state-space A matrix (symbolic or numeric)
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _parse_expr(expr_input: Any) -> sp.Expr:
    """Convert a string or number to a sympy expression."""
    if isinstance(expr_input, sp.Basic):
        return expr_input
    if isinstance(expr_input, (int, float)):
        return sp.Rational(expr_input) if isinstance(expr_input, int) else sp.nsimplify(expr_input)
    if isinstance(expr_input, str):
        return sp.sympify(expr_input, rational=True)
    return sp.sympify(expr_input)


def _classify_eigenvalues(eigenvalues: list[sp.Expr],
                          numeric: bool = False) -> dict[str, Any]:
    """Classify stability from a list of eigenvalues.

    Returns dict with:
        stability: str
        classification: str
        details: str
    """
    if numeric:
        # Numeric classification
        real_parts = [complex(ev).real for ev in eigenvalues]
        max_re = max(real_parts)
        min_re = min(real_parts)

        if max_re > 1e-10:
            stability = "Unstable"
            classification = "UNSTABLE"
            details = (
                f"At least one eigenvalue has positive real part "
                f"(max Re(lambda) = {max_re:.6g}). "
                f"The equilibrium is unstable -- perturbations grow exponentially."
            )
        elif max_re > -1e-10:
            # Some eigenvalues on imaginary axis
            on_axis = [ev for ev, rp in zip(eigenvalues, real_parts) if abs(rp) < 1e-10]
            # Check for repeated eigenvalues on imaginary axis
            from collections import Counter
            rounded = [round(complex(ev).imag, 8) for ev in on_axis]
            counts = Counter(rounded)
            has_repeated = any(v > 1 for v in counts.values())
            if has_repeated:
                stability = "Unstable"
                classification = "UNSTABLE"
                details = (
                    "Repeated eigenvalues on imaginary axis -- "
                    "secular (algebraic) growth."
                )
            else:
                stability = "Marginally stable"
                classification = "STABLE"
                details = (
                    "All eigenvalues have non-positive real part with "
                    "purely imaginary eigenvalues being simple. "
                    "The equilibrium is marginally stable (bounded oscillations)."
                )
        else:
            stability = "Asymptotically stable"
            classification = "STABLE"
            details = (
                f"All eigenvalues have strictly negative real parts "
                f"(max Re(lambda) = {max_re:.6g}). "
                f"Perturbations decay exponentially to zero."
            )
    else:
        # Symbolic -- try to determine sign
        stability = "Symbolic (see conditions)"
        classification = "CONDITIONAL"
        details = "Stability depends on parameter values."

    return {
        "stability": stability,
        "classification": classification,
        "details": details,
    }


def _sdof_stability(mass_coeff: sp.Expr,
                    damping_coeff: sp.Expr,
                    stiffness_coeff: sp.Expr,
                    param_values: dict | None = None) -> dict[str, Any]:
    """Analyse stability of SDOF: m*q'' + c*q' + k*q = 0.

    Returns dict with eigenvalues, stability classification, conditions, etc.
    """
    m = mass_coeff
    c = damping_coeff
    k = stiffness_coeff

    # State matrix: [[0, 1], [-k/m, -c/m]]
    lam = sp.Symbol("lambda")
    char_eq = lam**2 + (c / m) * lam + (k / m)

    eigenvalues_sym = sp.solve(char_eq, lam)

    result: dict[str, Any] = {
        "char_eq": char_eq,
        "eigenvalues_symbolic": eigenvalues_sym,
    }

    # Determine stability conditions symbolically
    # For m*q'' + c*q' + k*q = 0:
    #   effective_stiffness = k/m, effective_damping = c/m
    eff_k = sp.simplify(k / m)
    eff_c = sp.simplify(c / m)

    conditions = []
    if param_values is not None:
        # Numeric evaluation
        subs = {sp.Symbol(name): val for name, val in param_values.items()}
        k_num = float(eff_k.subs(subs))
        c_num = float(eff_c.subs(subs))
        eigs_num = [complex(ev.subs(subs)) for ev in eigenvalues_sym]
        result["eigenvalues_numeric"] = eigs_num
        result.update(_classify_eigenvalues(eigs_num, numeric=True))
    else:
        # Symbolic stability conditions
        conditions.append(f"For stability: k_eff = {eff_k} > 0 and c_eff = {eff_c} >= 0")
        conditions.append(f"  k_eff > 0 and c_eff > 0 => Asymptotically stable")
        conditions.append(f"  k_eff > 0 and c_eff = 0 => Marginally stable (undamped oscillation)")
        conditions.append(f"  k_eff < 0 => Unstable (divergent)")
        result["stability_conditions"] = conditions
        result["stability"] = "Symbolic (see conditions)"
        result["classification"] = "CONDITIONAL"
        result["details"] = "; ".join(conditions)

    result["effective_stiffness"] = eff_k
    result["effective_damping"] = eff_c
    return result


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------


class StabilitySolver(BaseSolver):
    """선형화된 시스템의 안정성 판별.

    Supports two input modes:
      1. "linearized_ode" -- SDOF coefficients (mass, damping, stiffness)
      2. "state_matrix"   -- general state-space A matrix
    """

    def solve(self, params: dict) -> SolverResult:
        """Perform stability analysis.

        Parameters:
            params: dict with keys depending on mode:

            Mode "linearized_ode":
                mass_coeff: str or sympy expr -- coefficient of q'' in linearized EOM
                damping_coeff: str or sympy expr -- coefficient of q'
                stiffness_coeff: str or sympy expr -- coefficient of q
                parameters: dict (optional) -- numeric values for symbols

            Mode "state_matrix":
                A: list of lists or sympy Matrix -- state matrix
                parameters: dict (optional) -- numeric values for symbols
        """
        mode = params.get("mode", "linearized_ode")

        if mode == "linearized_ode":
            return self._solve_linearized_ode(params)
        elif mode == "state_matrix":
            return self._solve_state_matrix(params)
        else:
            raise ValueError(f"Unknown stability mode: {mode}")

    def _solve_linearized_ode(self, params: dict) -> SolverResult:
        """Stability from SDOF linearized ODE coefficients."""
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        m_expr = _parse_expr(params.get("mass_coeff", "1"))
        c_expr = _parse_expr(params.get("damping_coeff", "0"))
        k_expr = _parse_expr(params.get("stiffness_coeff", "0"))
        param_values = params.get("parameters", None)

        # Step 1: System characterization
        is_damped = not sp.simplify(c_expr).is_zero if sp.simplify(c_expr).is_zero is not None else True
        sys_type = "Damped SDOF" if is_damped else "Undamped SDOF"
        steps.append((
            "System characterization",
            f"Type: {sys_type}\n"
            f"  EOM: ({m_expr})*q'' + ({c_expr})*q' + ({k_expr})*q = 0"
        ))

        # Step 2: Compute effective coefficients
        eff_k = sp.simplify(k_expr / m_expr)
        eff_c = sp.simplify(c_expr / m_expr)
        steps.append((
            "Effective coefficients",
            f"omega_n^2 = k_eff/m_eff = {eff_k}\n"
            f"2*zeta*omega_n = c_eff/m_eff = {eff_c}"
        ))

        # Step 3: Eigenvalue computation
        analysis = _sdof_stability(m_expr, c_expr, k_expr, param_values)
        eigs_sym = analysis["eigenvalues_symbolic"]
        eigs_str = ", ".join(str(e) for e in eigs_sym)
        steps.append((
            "Eigenvalue computation (symbolic)",
            f"Characteristic equation: lambda^2 + ({eff_c})*lambda + ({eff_k}) = 0\n"
            f"Eigenvalues: {eigs_str}"
        ))

        if param_values is not None:
            eigs_num = analysis.get("eigenvalues_numeric", [])
            eigs_num_str = ", ".join(f"{e}" for e in eigs_num)
            steps.append((
                "Eigenvalue computation (numeric)",
                f"Eigenvalues: {eigs_num_str}"
            ))

        # Step 4: Stability classification
        stability = analysis["stability"]
        classification = analysis["classification"]
        details = analysis["details"]
        steps.append((
            "Stability classification",
            f"Result: {stability}\n{details}"
        ))

        # Step 5: Physical interpretation
        if param_values is not None:
            subs = {sp.Symbol(name): val for name, val in param_values.items()}
            k_val = float(eff_k.subs(subs))
            if k_val > 0:
                omega_n = float(sp.sqrt(eff_k).subs(subs))
                interp = (
                    f"Effective stiffness k_eff = {k_val:.6g} > 0.\n"
                    f"Natural frequency omega_n = {omega_n:.6g} rad/s.\n"
                    f"System oscillates about equilibrium."
                )
            elif k_val < 0:
                growth_rate = float(sp.sqrt(-eff_k).subs(subs))
                interp = (
                    f"Effective stiffness k_eff = {k_val:.6g} < 0.\n"
                    f"Growth rate = {growth_rate:.6g} (1/s).\n"
                    f"System diverges exponentially from equilibrium."
                )
            else:
                interp = "Effective stiffness = 0: neutral equilibrium."
            steps.append(("Physical interpretation", interp))
        else:
            # Symbolic stability condition
            cond_str = (
                f"Stability condition: {eff_k} > 0\n"
                f"If additionally {eff_c} > 0: asymptotically stable.\n"
                f"If {eff_c} = 0 and {eff_k} > 0: marginally stable."
            )
            steps.append(("Stability condition", cond_str))
            final_answer["stability_condition"] = f"{eff_k} > 0"

        final_answer["stability"] = stability
        final_answer["classification"] = classification
        final_answer["eigenvalues_symbolic"] = [str(e) for e in eigs_sym]
        final_answer["effective_stiffness"] = str(eff_k)
        final_answer["effective_damping"] = str(eff_c)
        if param_values is not None and "eigenvalues_numeric" in analysis:
            final_answer["eigenvalues_numeric"] = [
                complex(e) for e in analysis["eigenvalues_numeric"]
            ]

        # Sanity check
        sanity_parts = []
        sanity_parts.append(
            f"Stability verdict: {stability} ({classification})"
        )
        if param_values is not None:
            subs = {sp.Symbol(name): val for name, val in param_values.items()}
            k_val = float(eff_k.subs(subs))
            sanity_parts.append(f"k_eff = {k_val:.6g} ({'> 0 => restoring' if k_val > 0 else '< 0 => divergent'})")

        return SolverResult(
            problem_type="Stability Analysis (Linearized ODE)",
            given={
                "mass_coeff": str(m_expr),
                "damping_coeff": str(c_expr),
                "stiffness_coeff": str(k_expr),
                "parameters": param_values,
            },
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    def _solve_state_matrix(self, params: dict) -> SolverResult:
        """Stability from state-space A matrix."""
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        A_input = params.get("A")
        if A_input is None:
            raise ValueError("'A' (state matrix) is required for state_matrix mode.")

        param_values = params.get("parameters", None)

        # Build sympy Matrix
        if isinstance(A_input, sp.Matrix):
            A = A_input
        else:
            A = sp.Matrix(A_input)

        n = A.shape[0]
        steps.append((
            "System characterization",
            f"State-space system: dx/dt = A*x\n"
            f"A matrix ({n}x{n}):\n{sp.sstr(A)}"
        ))

        # Eigenvalues
        lam = sp.Symbol("lambda")
        char_poly = A.charpoly(lam)
        char_eq = char_poly.as_expr()
        eigenvalues = A.eigenvals()  # dict: eigenvalue -> multiplicity

        eig_list = []
        for ev, mult in eigenvalues.items():
            for _ in range(mult):
                eig_list.append(ev)

        eigs_str = ", ".join(str(e) for e in eig_list)
        steps.append((
            "Characteristic polynomial",
            f"det(A - lambda*I) = {char_eq}\n"
        ))
        steps.append((
            "Eigenvalue computation",
            f"Eigenvalues (with multiplicity): {eigs_str}"
        ))

        # Numeric evaluation if parameters given
        if param_values is not None:
            subs = {sp.Symbol(name): val for name, val in param_values.items()}
            eigs_num = [complex(ev.subs(subs)) for ev in eig_list]
            eigs_num_str = ", ".join(f"{e}" for e in eigs_num)
            steps.append((
                "Eigenvalues (numeric)",
                eigs_num_str,
            ))
            cls = _classify_eigenvalues(eigs_num, numeric=True)
        else:
            # Try numeric if all entries are numbers
            try:
                eigs_num = [complex(ev) for ev in eig_list]
                cls = _classify_eigenvalues(eigs_num, numeric=True)
            except (TypeError, ValueError):
                cls = _classify_eigenvalues(eig_list, numeric=False)

        steps.append((
            "Stability classification",
            f"Result: {cls['stability']}\n{cls['details']}"
        ))

        final_answer["stability"] = cls["stability"]
        final_answer["classification"] = cls["classification"]
        final_answer["eigenvalues"] = [str(e) for e in eig_list]
        final_answer["characteristic_polynomial"] = str(char_eq)

        return SolverResult(
            problem_type="Stability Analysis (State Matrix)",
            given={"A": str(A), "parameters": param_values},
            steps=steps,
            final_answer=final_answer,
            sanity_check=f"Stability verdict: {cls['stability']} ({cls['classification']})",
        )

    def get_input_template(self) -> dict:
        """필요한 입력 파라미터 설명을 반환한다."""
        return {
            "mode": {
                "type": "str",
                "required": False,
                "default": "linearized_ode",
                "description": "Analysis mode: 'linearized_ode' or 'state_matrix'",
            },
            "mass_coeff": {
                "type": "str or sympy expr",
                "required": "for linearized_ode mode",
                "description": "Coefficient of q'' in the linearized EOM",
            },
            "damping_coeff": {
                "type": "str or sympy expr",
                "required": "for linearized_ode mode",
                "description": "Coefficient of q' in the linearized EOM",
            },
            "stiffness_coeff": {
                "type": "str or sympy expr",
                "required": "for linearized_ode mode",
                "description": "Coefficient of q in the linearized EOM",
            },
            "A": {
                "type": "list of lists or sympy Matrix",
                "required": "for state_matrix mode",
                "description": "State matrix A for dx/dt = Ax",
            },
            "parameters": {
                "type": "dict",
                "required": False,
                "description": "Numeric parameter values for substitution, e.g. {'m': 1, 'g': 9.81}",
            },
        }
