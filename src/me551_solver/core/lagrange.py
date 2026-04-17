"""Lagrange 방정식 -> 평형점 -> 선형화 -> 안정성 분석 솔버.

Symbolically derives equations of motion using Lagrange's equation.
Supports 1-DOF (legacy) and N-DOF (v2) systems.
"""

from __future__ import annotations

from typing import Any

import sympy as sp
from sympy import (
    Derivative, Function, Rational, Symbol, cos, diff, expand,
    series, simplify, sin, solve, sqrt, symbols, trigsimp,
)

from .base import BaseSolver, SolverResult
from .stability import StabilitySolver


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _parse_expr(expr_input: Any, local_dict: dict | None = None) -> sp.Expr:
    """Parse a string or numeric value into a sympy expression."""
    if isinstance(expr_input, sp.Basic):
        return expr_input
    if expr_input is None:
        return sp.S.Zero
    if isinstance(expr_input, (int, float)):
        return sp.nsimplify(expr_input)
    if isinstance(expr_input, str):
        return sp.sympify(expr_input, locals=local_dict, rational=True)
    return sp.sympify(expr_input)


def _build_local_dict(symbol_names: list[str],
                      coord: str,
                      coord_dot: str) -> dict[str, Any]:
    """Build a local dictionary for sympify parsing.

    Creates Symbol objects for all parameter names, the coordinate,
    and the coordinate velocity.
    """
    local = {}
    # Standard sympy functions and helpers
    local["Rational"] = sp.Rational
    local["cos"] = sp.cos
    local["sin"] = sp.sin
    local["sqrt"] = sp.sqrt
    local["pi"] = sp.pi
    local["tan"] = sp.tan
    local["atan"] = sp.atan
    local["atan2"] = sp.atan2
    local["exp"] = sp.exp
    local["log"] = sp.log
    local["Abs"] = sp.Abs

    # Create symbols
    for name in symbol_names:
        if name not in local:
            local[name] = sp.Symbol(name, real=True)

    # Ensure coord and coord_dot are symbols
    if coord not in local:
        local[coord] = sp.Symbol(coord, real=True)
    if coord_dot not in local:
        local[coord_dot] = sp.Symbol(coord_dot, real=True)

    return local


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------


class LagrangeSolver(BaseSolver):
    """Lagrange 방정식 유도, 평형점 계산, 선형화, 안정성 판별.

    Accepts two input formats:
      1. Detailed format with "coordinates", "kinetic_energy", etc.
      2. Simpler format with "T_expr", "V_expr", "coord", "coord_dot", etc.
    """

    def solve(self, params: dict) -> SolverResult:
        """Derive EOM via Lagrange's equation and perform stability analysis.

        Automatically detects 1-DOF vs N-DOF from input format.
        For N-DOF, use ``coords`` as a list of (q, qdot) pairs.
        """
        # Check for multi-DOF
        if "coords" in params and len(params["coords"]) > 1:
            return self._solve_ndof(params)

        # Legacy 1-DOF path
        parsed = self._parse_input(params)

        coord_name = parsed["coord"]
        coord_dot_name = parsed["coord_dot"]
        T_str = parsed["T_expr"]
        V_str = parsed["V_expr"]
        Qnc_str = parsed.get("Qnc_expr", "0")
        symbol_names = parsed["symbols"]
        find = parsed.get("find", ["eom", "equilibrium", "linearization", "stability"])

        # Build symbol dictionary
        local = _build_local_dict(symbol_names, coord_name, coord_dot_name)
        q_sym = local[coord_name]      # e.g., Symbol('r')
        qdot_sym = local[coord_dot_name]  # e.g., Symbol('rdot')

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # ------------------------------------------------------------------
        # Step 1: Define generalized coordinate
        # ------------------------------------------------------------------
        steps.append((
            "Define generalized coordinate",
            f"q = {coord_name}, q_dot = {coord_dot_name}"
        ))

        # ------------------------------------------------------------------
        # Step 2: Parse kinetic energy T
        # ------------------------------------------------------------------
        T = _parse_expr(T_str, local)
        steps.append((
            "Kinetic energy T",
            f"T = {T}"
        ))
        final_answer["T"] = str(T)

        # ------------------------------------------------------------------
        # Step 3: Parse potential energy V
        # ------------------------------------------------------------------
        V = _parse_expr(V_str, local)
        steps.append((
            "Potential energy V",
            f"V = {V}"
        ))
        final_answer["V"] = str(V)

        # ------------------------------------------------------------------
        # Step 4: Compute dT/d(qdot)
        # ------------------------------------------------------------------
        dT_dqdot = sp.diff(T, qdot_sym)
        dT_dqdot = sp.simplify(dT_dqdot)
        steps.append((
            f"Compute dT/d({coord_dot_name})",
            f"dT/d({coord_dot_name}) = {dT_dqdot}"
        ))
        final_answer["dT_dqdot"] = str(dT_dqdot)

        # ------------------------------------------------------------------
        # Step 5: Compute d/dt(dT/dqdot)
        # For static expressions (no explicit t-dependence in T written in
        # terms of q and qdot), d/dt(dT/dqdot) is obtained by:
        #   replacing qdot -> qddot (coefficient of qdot in dT/dqdot)
        #   plus chain-rule terms from q -> qdot
        # We use the shortcut: dT/dqdot is linear in qdot for standard
        # kinetic energies, so d/dt(dT/dqdot) = coeff_qdot * qddot
        #   + (dT/dqdot w.r.t. q) * qdot
        # ------------------------------------------------------------------
        qddot_sym = sp.Symbol(coord_name + "ddot", real=True)

        # d/dt of dT/dqdot: replace qdot->qddot via chain rule
        # d/dt(f(q, qdot)) = df/dq * qdot + df/dqdot * qddot
        ddt_dT_dqdot = (sp.diff(dT_dqdot, q_sym) * qdot_sym
                        + sp.diff(dT_dqdot, qdot_sym) * qddot_sym)
        ddt_dT_dqdot = sp.simplify(ddt_dT_dqdot)
        steps.append((
            f"Compute d/dt(dT/d({coord_dot_name}))",
            f"d/dt(dT/d({coord_dot_name})) = {ddt_dT_dqdot}"
        ))
        final_answer["ddt_dT_dqdot"] = str(ddt_dT_dqdot)

        # ------------------------------------------------------------------
        # Step 6: Compute dT/dq
        # ------------------------------------------------------------------
        dT_dq = sp.diff(T, q_sym)
        dT_dq = sp.simplify(dT_dq)
        steps.append((
            f"Compute dT/d{coord_name}",
            f"dT/d{coord_name} = {dT_dq}"
        ))
        final_answer["dT_dq"] = str(dT_dq)

        # ------------------------------------------------------------------
        # Step 7: Compute dV/dq
        # ------------------------------------------------------------------
        dV_dq = sp.diff(V, q_sym)
        dV_dq = sp.simplify(dV_dq)
        steps.append((
            f"Compute dV/d{coord_name}",
            f"dV/d{coord_name} = {dV_dq}"
        ))
        final_answer["dV_dq"] = str(dV_dq)

        # ------------------------------------------------------------------
        # Step 8: Lagrange equation:
        #   d/dt(dT/dqdot) - dT/dq + dV/dq = Q_nc
        # ------------------------------------------------------------------
        Qnc = _parse_expr(Qnc_str, local)
        eom = sp.simplify(ddt_dT_dqdot - dT_dq + dV_dq - Qnc)
        eom = sp.expand(eom)

        steps.append((
            f"Lagrange equation: d/dt(dT/d{coord_dot_name}) - dT/d{coord_name} + dV/d{coord_name} = Q_nc",
            f"EOM: {eom} = 0"
        ))
        final_answer["eom"] = str(eom)

        # Extract mass coefficient (coefficient of qddot)
        mass_coeff = sp.simplify(eom.coeff(qddot_sym))
        if mass_coeff == 0:
            mass_coeff = sp.S.One
        final_answer["mass_coeff"] = str(mass_coeff)

        # ------------------------------------------------------------------
        # Step 9: Equilibrium points (set qdot = 0, qddot = 0)
        # ------------------------------------------------------------------
        if "equilibrium" in find:
            eom_static = eom.subs(qdot_sym, 0).subs(qddot_sym, 0)
            eom_static = sp.simplify(eom_static)

            # Try to solve for q
            eq_solutions = sp.solve(eom_static, q_sym)

            steps.append((
                "Equilibrium points",
                f"Set {coord_dot_name} = 0, {coord_name}ddot = 0:\n"
                f"  {eom_static} = 0\n"
                f"  Solutions: {coord_name} = {eq_solutions}"
            ))
            final_answer["equilibrium_equation"] = str(eom_static)
            final_answer["equilibrium_points"] = [str(s) for s in eq_solutions]

        # ------------------------------------------------------------------
        # Step 10 & 11: Linearization around equilibrium
        # ------------------------------------------------------------------
        if "linearization" in find and "equilibrium" in find:
            delta = sp.Symbol("delta_" + coord_name, real=True)
            delta_dot = sp.Symbol("delta_" + coord_dot_name, real=True)
            delta_ddot = sp.Symbol("delta_" + coord_name + "ddot", real=True)

            linearized_results = []

            for eq_pt in eq_solutions:
                steps.append((
                    f"Introduce perturbation about {coord_name}_eq = {eq_pt}",
                    f"{coord_name} = {eq_pt} + delta_{coord_name}\n"
                    f"{coord_dot_name} = delta_{coord_dot_name}\n"
                    f"{coord_name}ddot = delta_{coord_name}ddot"
                ))

                # Substitute q = q_eq + delta into the EOM
                eom_perturbed = eom.subs(q_sym, eq_pt + delta)
                eom_perturbed = eom_perturbed.subs(qdot_sym, delta_dot)
                eom_perturbed = eom_perturbed.subs(qddot_sym, delta_ddot)

                # Taylor expand: keep terms up to first order in delta
                # Use series expansion around delta=0
                eom_lin = sp.S.Zero
                # Constant term (should be zero at equilibrium)
                const_term = eom_perturbed.subs(delta, 0).subs(delta_dot, 0).subs(delta_ddot, 0)

                # Coefficient of delta_ddot
                coeff_ddot = sp.simplify(eom_perturbed.coeff(delta_ddot))

                # Coefficient of delta_dot
                coeff_dot = sp.simplify(eom_perturbed.coeff(delta_dot))

                # Coefficient of delta (first derivative of EOM w.r.t. delta at delta=0)
                # Remove ddot and dot terms first, then differentiate
                eom_no_accel = eom_perturbed - coeff_ddot * delta_ddot - coeff_dot * delta_dot
                coeff_delta = sp.simplify(sp.diff(eom_no_accel, delta).subs(delta, 0))

                eom_lin = coeff_ddot * delta_ddot + coeff_dot * delta_dot + coeff_delta * delta

                # Simplify
                eom_lin = sp.simplify(eom_lin)
                coeff_ddot = sp.simplify(coeff_ddot.subs(delta, 0).subs(delta_dot, 0))
                coeff_dot = sp.simplify(coeff_dot.subs(delta, 0).subs(delta_dot, 0))
                coeff_delta = sp.simplify(coeff_delta)

                steps.append((
                    f"Linearized EOM about {coord_name}_eq = {eq_pt}",
                    f"({coeff_ddot})*delta_{coord_name}ddot + "
                    f"({coeff_dot})*delta_{coord_dot_name} + "
                    f"({coeff_delta})*delta_{coord_name} = 0"
                ))

                steps.append((
                    f"Effective stiffness about {coord_name}_eq = {eq_pt}",
                    f"k_eff = {coeff_delta}\n"
                    f"m_eff = {coeff_ddot}\n"
                    f"c_eff = {coeff_dot}"
                ))

                lin_info = {
                    "equilibrium": str(eq_pt),
                    "mass_coeff": str(coeff_ddot),
                    "damping_coeff": str(coeff_dot),
                    "stiffness_coeff": str(coeff_delta),
                }
                linearized_results.append(lin_info)

            final_answer["linearized"] = linearized_results

        # ------------------------------------------------------------------
        # Step 12 & 13: Stability classification
        # ------------------------------------------------------------------
        if "stability" in find and "linearization" in find and "equilibrium" in find:
            stability_solver = StabilitySolver()
            param_values = parsed.get("parameters", None)
            stability_results = []

            for lin_info in linearized_results:
                stab_params: dict[str, Any] = {
                    "mode": "linearized_ode",
                    "mass_coeff": lin_info["mass_coeff"],
                    "damping_coeff": lin_info["damping_coeff"],
                    "stiffness_coeff": lin_info["stiffness_coeff"],
                }
                if param_values:
                    stab_params["parameters"] = param_values

                stab_result = stability_solver.solve(stab_params)

                eq_pt = lin_info["equilibrium"]
                stability_str = stab_result.final_answer.get("stability", "Unknown")
                classification = stab_result.final_answer.get("classification", "Unknown")

                steps.append((
                    f"Stability classification for {coord_name}_eq = {eq_pt}",
                    f"Stability: {stability_str} ({classification})\n"
                    + "\n".join(f"  {desc}: {val}" for desc, val in stab_result.steps)
                ))

                stab_info = {
                    "equilibrium": eq_pt,
                    "stability": stability_str,
                    "classification": classification,
                    "eigenvalues": stab_result.final_answer.get("eigenvalues_symbolic", []),
                }
                if "stability_condition" in stab_result.final_answer:
                    stab_info["stability_condition"] = stab_result.final_answer["stability_condition"]
                if "eigenvalues_numeric" in stab_result.final_answer:
                    stab_info["eigenvalues_numeric"] = [
                        str(e) for e in stab_result.final_answer["eigenvalues_numeric"]
                    ]
                stability_results.append(stab_info)

            final_answer["stability"] = stability_results

            # Set top-level stability for single-equilibrium convenience
            if len(stability_results) == 1:
                final_answer["stability_verdict"] = stability_results[0]["stability"]
                final_answer["classification"] = stability_results[0]["classification"]

        # ------------------------------------------------------------------
        # Sanity check
        # ------------------------------------------------------------------
        sanity_parts = []
        # Verify EOM is second-order in acceleration
        if mass_coeff != 0:
            sanity_parts.append(f"Mass coefficient = {mass_coeff} (non-zero => valid 2nd-order ODE)")
        if "equilibrium" in find:
            sanity_parts.append(f"Number of equilibrium points found: {len(eq_solutions)}")
        if "stability" in find and "stability" in final_answer:
            for sr in final_answer["stability"]:
                sanity_parts.append(
                    f"Equilibrium {coord_name}={sr['equilibrium']}: {sr['stability']}"
                )

        return SolverResult(
            problem_type="Lagrange Equation Analysis",
            given={
                "T": T_str,
                "V": V_str,
                "coordinate": coord_name,
                "symbols": symbol_names,
            },
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    # ------------------------------------------------------------------
    # N-DOF Lagrange solver (v2)
    # ------------------------------------------------------------------

    def _solve_ndof(self, params: dict) -> SolverResult:
        """Derive EOM for an N-DOF system via Lagrange's equations.

        Parameters
        ----------
        params : dict
            - ``coords``: list of (q_name, qdot_name) pairs
            - ``T_expr``: kinetic energy string
            - ``V_expr``: potential energy string
            - ``parameters``: list of symbol names
            - ``Qnc_exprs``: dict {q_name: expr_str} for non-conservative forces
            - ``parameter_values``: optional numeric values
            - ``find``: list of stages
        """
        coords = params["coords"]  # e.g. [("q1","q1dot"), ("q2","q2dot")]
        n_dof = len(coords)
        T_str = params.get("T_expr", "0")
        V_str = params.get("V_expr", "0")
        Qnc_strs = params.get("Qnc_exprs", {})
        symbol_names = params.get("parameters", [])
        find = params.get("find", ["eom", "equilibrium", "linearization", "stability"])
        param_values = params.get("parameter_values", None)

        # Build symbol dictionary
        local: dict[str, Any] = {}
        local["Rational"] = sp.Rational
        for fn_name in ("cos", "sin", "sqrt", "tan", "atan", "atan2", "exp", "log", "Abs"):
            local[fn_name] = getattr(sp, fn_name)
        local["pi"] = sp.pi
        for name in symbol_names:
            if name not in local:
                local[name] = sp.Symbol(name, real=True)

        q_syms: list[sp.Symbol] = []
        qdot_syms: list[sp.Symbol] = []
        qddot_syms: list[sp.Symbol] = []
        for q_name, qdot_name in coords:
            q_s = sp.Symbol(q_name, real=True)
            qd_s = sp.Symbol(qdot_name, real=True)
            qdd_s = sp.Symbol(q_name + "ddot", real=True)
            q_syms.append(q_s)
            qdot_syms.append(qd_s)
            qddot_syms.append(qdd_s)
            local[q_name] = q_s
            local[qdot_name] = qd_s
            local[q_name + "ddot"] = qdd_s

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # Step 1: coordinates
        coord_strs = [f"{q_name} / {qdot_name}" for q_name, qdot_name in coords]
        steps.append(("Generalized coordinates", f"{n_dof}-DOF: " + ", ".join(coord_strs)))

        # Step 2-3: Parse T, V
        T = _parse_expr(T_str, local)
        V = _parse_expr(V_str, local)
        steps.append(("Kinetic energy T", f"T = {T}"))
        steps.append(("Potential energy V", f"V = {V}"))
        final_answer["T"] = str(T)
        final_answer["V"] = str(V)

        # Step 4: Derive EOM for each coordinate
        eom_list: list[sp.Expr] = []
        for i in range(n_dof):
            q_i = q_syms[i]
            qd_i = qdot_syms[i]
            qdd_i = qddot_syms[i]

            # dT/d(qdot_i)
            dT_dqdi = sp.simplify(sp.diff(T, qd_i))

            # d/dt(dT/dqdot_i) via chain rule
            ddt = sp.S.Zero
            for j in range(n_dof):
                ddt += sp.diff(dT_dqdi, q_syms[j]) * qdot_syms[j]
                ddt += sp.diff(dT_dqdi, qdot_syms[j]) * qddot_syms[j]
            ddt = sp.simplify(ddt)

            # dT/dq_i
            dT_dqi = sp.simplify(sp.diff(T, q_i))

            # dV/dq_i
            dV_dqi = sp.simplify(sp.diff(V, q_i))

            # Non-conservative force
            q_name = coords[i][0]
            Qnc_i = _parse_expr(Qnc_strs.get(q_name, "0"), local)

            eom_i = sp.simplify(sp.expand(ddt - dT_dqi + dV_dqi - Qnc_i))
            eom_list.append(eom_i)

            steps.append((
                f"EOM for {q_name}",
                f"d/dt(dT/d{coords[i][1]}) - dT/d{q_name} + dV/d{q_name} = Q_nc\n"
                f"=> {eom_i} = 0"
            ))

        final_answer["eom"] = [str(e) for e in eom_list]

        # Step 5: Equilibrium
        eq_solutions = []
        if "equilibrium" in find:
            # Set all velocities and accelerations to zero
            subs_zero = {}
            for qd_i in qdot_syms:
                subs_zero[qd_i] = 0
            for qdd_i in qddot_syms:
                subs_zero[qdd_i] = 0

            static_eoms = [sp.simplify(e.subs(subs_zero)) for e in eom_list]
            eq_solutions = sp.solve(static_eoms, q_syms, dict=True)

            if not eq_solutions:
                # Try individually
                eq_solutions = [{}]
                for i, eq in enumerate(static_eoms):
                    sol = sp.solve(eq, q_syms[i])
                    if sol:
                        eq_solutions[0][q_syms[i]] = sol[0]

            steps.append((
                "Equilibrium points",
                "\n".join(f"  {e} = 0" for e in static_eoms)
                + f"\n  Solutions: {eq_solutions}"
            ))
            final_answer["equilibrium_points"] = [
                {str(k): str(v) for k, v in sol.items()} for sol in eq_solutions
            ]

        # Step 6: Linearization via Jacobian
        linearized_results = []
        if "linearization" in find and eq_solutions:
            for eq_pt in eq_solutions:
                # Substitute equilibrium into EOM
                subs_eq = {}
                subs_eq.update({qd: 0 for qd in qdot_syms})
                subs_eq.update(eq_pt)

                # Extract M, C, K matrices via coefficient extraction
                M_mat = sp.zeros(n_dof)
                C_mat = sp.zeros(n_dof)
                K_mat = sp.zeros(n_dof)

                for i in range(n_dof):
                    eom_expanded_i = sp.expand(eom_list[i])
                    for j in range(n_dof):
                        # Mass: coefficient of qddot_j in eom_i
                        M_mat[i, j] = sp.simplify(
                            eom_expanded_i.coeff(qddot_syms[j]).subs(subs_eq)
                        )
                        # Damping: coefficient of qdot_j in eom_i
                        C_mat[i, j] = sp.simplify(
                            eom_expanded_i.coeff(qdot_syms[j]).subs(subs_eq)
                        )

                # Stiffness: Jacobian of static EOM w.r.t. q
                subs_zero_vel = {qd: 0 for qd in qdot_syms}
                subs_zero_vel.update({qdd: 0 for qdd in qddot_syms})
                for i in range(n_dof):
                    eom_static_i = eom_list[i].subs(subs_zero_vel)
                    for j in range(n_dof):
                        K_mat[i, j] = sp.simplify(
                            sp.diff(eom_static_i, q_syms[j]).subs(eq_pt)
                        )

                eq_str = ", ".join(f"{k}={v}" for k, v in eq_pt.items())
                steps.append((
                    f"Linearization about ({eq_str})",
                    f"M = {M_mat}\n"
                    f"C = {C_mat}\n"
                    f"K = {K_mat}"
                ))

                lin_info = {
                    "equilibrium": {str(k): str(v) for k, v in eq_pt.items()},
                    "M": str(M_mat.tolist()),
                    "C": str(C_mat.tolist()),
                    "K": str(K_mat.tolist()),
                }

                # Numeric evaluation if parameter values given
                if param_values:
                    subs_num = {sp.Symbol(k, real=True): v for k, v in param_values.items()}
                    M_num = M_mat.subs(subs_num)
                    C_num = C_mat.subs(subs_num)
                    K_num = K_mat.subs(subs_num)
                    lin_info["M_numeric"] = str(M_num.tolist())
                    lin_info["C_numeric"] = str(C_num.tolist())
                    lin_info["K_numeric"] = str(K_num.tolist())
                    steps.append((
                        f"Numeric matrices at ({eq_str})",
                        f"M = {M_num}\nC = {C_num}\nK = {K_num}"
                    ))

                linearized_results.append(lin_info)

            final_answer["linearized"] = linearized_results

        # Step 7: Stability
        if "stability" in find and linearized_results:
            stability_results = []
            for idx_lin, lin_info in enumerate(linearized_results):
                try:
                    import numpy as np

                    subs_num = {}
                    if param_values:
                        subs_num = {sp.Symbol(k, real=True): v
                                    for k, v in param_values.items()}

                    # Use the sympy matrices we already computed, not string repr
                    eq_idx = idx_lin
                    eq_pt = eq_solutions[eq_idx] if eq_idx < len(eq_solutions) else {}
                    subs_eq_num = {}
                    subs_eq_num.update({qd: 0 for qd in qdot_syms})
                    subs_eq_num.update(eq_pt)

                    # Recompute M, K numerically
                    M_recomp = sp.zeros(n_dof)
                    K_recomp = sp.zeros(n_dof)
                    subs_zero_vel = {qd: 0 for qd in qdot_syms}
                    subs_zero_vel.update({qdd: 0 for qdd in qddot_syms})
                    for i in range(n_dof):
                        for j in range(n_dof):
                            M_recomp[i, j] = eom_list[i].coeff(qddot_syms[j]).subs(subs_eq_num)
                        eom_static_i = eom_list[i].subs(subs_zero_vel)
                        for j in range(n_dof):
                            K_recomp[i, j] = sp.diff(eom_static_i, q_syms[j]).subs(eq_pt)

                    M_num = np.array(M_recomp.subs(subs_num).tolist(), dtype=float)
                    K_num = np.array(K_recomp.subs(subs_num).tolist(), dtype=float)

                    # Eigenvalue problem: K u = λ M u
                    from scipy.linalg import eig
                    eigenvalues, _ = eig(K_num, M_num)
                    eigenvalues = eigenvalues.real

                    # Stability from K eigenvalues (for conservative system)
                    all_positive = all(ev > -1e-10 for ev in eigenvalues)
                    has_negative = any(ev < -1e-10 for ev in eigenvalues)

                    if has_negative:
                        verdict = "UNSTABLE"
                    elif all_positive:
                        verdict = "STABLE (marginally)"
                    else:
                        verdict = "CONDITIONAL"

                    eq_str = ", ".join(
                        f"{k}={v}" for k, v in lin_info["equilibrium"].items()
                    )
                    steps.append((
                        f"Stability at ({eq_str})",
                        f"K eigenvalues: {eigenvalues.tolist()}\n"
                        f"Verdict: {verdict}"
                    ))

                    stability_results.append({
                        "equilibrium": lin_info["equilibrium"],
                        "stability": verdict,
                        "K_eigenvalues": eigenvalues.tolist(),
                    })
                except Exception as exc:
                    # Symbolic stability (no numeric values)
                    stability_results.append({
                        "equilibrium": lin_info["equilibrium"],
                        "stability": "SYMBOLIC (numeric evaluation needed)",
                        "note": str(exc),
                    })

            final_answer["stability"] = stability_results

        # Sanity check
        sanity_parts = [f"{n_dof}-DOF Lagrange analysis"]
        if eq_solutions:
            sanity_parts.append(f"Equilibrium points found: {len(eq_solutions)}")
        if "stability" in final_answer:
            for sr in final_answer["stability"]:
                sanity_parts.append(
                    f"  {sr['equilibrium']}: {sr['stability']}"
                )

        return SolverResult(
            problem_type="Lagrange Equation Analysis (N-DOF)",
            given={
                "T": T_str,
                "V": V_str,
                "coordinates": [c[0] for c in coords],
                "n_dof": n_dof,
                "symbols": symbol_names,
            },
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    def _parse_input(self, params: dict) -> dict:
        """Normalize both input formats to a common internal representation."""
        result: dict[str, Any] = {}

        if "T_expr" in params:
            # Simpler format
            result["coord"] = params.get("coord", "q")
            result["coord_dot"] = params.get("coord_dot", result["coord"] + "dot")
            result["T_expr"] = params["T_expr"]
            result["V_expr"] = params.get("V_expr", "0")
            result["Qnc_expr"] = params.get("Qnc_expr", "0")
            result["symbols"] = params.get("parameters", [])
            result["find"] = params.get("find", ["eom", "equilibrium", "linearization", "stability"])
            result["parameters"] = params.get("parameter_values", None)
        elif "kinetic_energy" in params:
            # Detailed format
            coords = params.get("coordinates", {"q": "q"})
            # For 1-DOF, take first coordinate
            coord_label = list(coords.values())[0] if coords else "q"
            result["coord"] = coord_label
            result["coord_dot"] = coord_label + "dot"
            result["T_expr"] = params["kinetic_energy"]
            result["V_expr"] = params.get("potential_energy", "0")
            # Non-conservative forces
            nc_forces = params.get("nonconservative_forces", {})
            result["Qnc_expr"] = nc_forces.get(coord_label, "0") if nc_forces else "0"
            result["symbols"] = params.get("symbols", [])
            result["find"] = params.get("find", ["eom", "equilibrium", "linearization", "stability"])
            result["parameters"] = params.get("parameter_values", None)
        else:
            raise ValueError(
                "Input must contain either 'T_expr' (simple format) or "
                "'kinetic_energy' (detailed format)."
            )

        return result

    def get_input_template(self) -> dict:
        """필요한 입력 파라미터 설명을 반환한다."""
        return {
            "T_expr": {
                "type": "str",
                "required": True,
                "description": (
                    "Kinetic energy T as a sympy-parseable string expression. "
                    "Use Rational(1,2) for fractions. "
                    "Example: 'Rational(1,2)*m*(rdot**2 + (Omega*r)**2)'"
                ),
            },
            "V_expr": {
                "type": "str",
                "required": False,
                "default": "0",
                "description": (
                    "Potential energy V as a sympy-parseable string. "
                    "Example: 'm*g*r'"
                ),
            },
            "coord": {
                "type": "str",
                "required": False,
                "default": "q",
                "description": "Name of the generalized coordinate (e.g. 'r', 'theta').",
            },
            "coord_dot": {
                "type": "str",
                "required": False,
                "default": "qdot",
                "description": "Name of the generalized velocity (e.g. 'rdot', 'thetadot').",
            },
            "parameters": {
                "type": "list[str]",
                "required": False,
                "description": (
                    "List of symbolic parameter names. "
                    "Example: ['m', 'g', 'Omega']"
                ),
            },
            "parameter_values": {
                "type": "dict",
                "required": False,
                "description": (
                    "Numeric values for parameters (optional). "
                    "Example: {'m': 1, 'g': 9.81, 'Omega': 2}"
                ),
            },
            "find": {
                "type": "list[str]",
                "required": False,
                "default": ["eom", "equilibrium", "linearization", "stability"],
                "description": (
                    "What to compute: 'eom', 'equilibrium', 'linearization', 'stability'."
                ),
            },
        }
