"""2-DOF 체인(직렬 스프링-질량) 문제 템플릿.

Automatically assembles the mass matrix M and stiffness matrix K from
physical parameters and delegates the eigenvalue / modal analysis to
``core.modal.ModalSolver`` once that is implemented.  In the mean time
the assembly and full step-by-step derivation are provided here using
NumPy/SymPy so that the template is immediately useful.

Supported topologies
--------------------
"chain" (default)
    wall ─ k1 ─ [m1] ─ k2 ─ [m2] ─ k3 ─ wall

    M = diag(m1, m2)
    K = [[k1+k2,  -k2  ],
         [-k2,    k2+k3]]

    k3 = 0  →  free right end
    k1 = 0  →  free left end (unusual but accepted)

"grounded_both_ends"
    Alias for "chain" with k1 > 0 and k3 > 0 (explicit).

"free_free"
    [m1] ─ k2 ─ [m2]   (no wall springs, k1=k3=0)

    M = diag(m1, m2)
    K = [[k2, -k2],
         [-k2, k2]]

    This configuration has a rigid-body mode (zero natural frequency).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import sympy as sp

from ..core.base import BaseSolver, SolverResult


class TwoDOFChainSolver(BaseSolver):
    """직렬 연결 2-DOF 스프링-질량 시스템: 행렬 조립 → 모달 해석.

    Supports symbolic (None) and numeric parameters.
    """

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self, params: dict) -> SolverResult:
        """Solve the 2-DOF chain system.

        Parameters
        ----------
        params : dict
            masses   : list[float | None]  – [m1, m2]
            springs  : list[float | None]  – [k1, k2] or [k1, k2, k3]
                        k1 = wall-to-m1, k2 = m1-to-m2, k3 = m2-to-wall
                        (k3 defaults to 0 if only two springs given)
            topology : str  – "chain", "grounded_both_ends", or "free_free"
            initial_conditions : dict | None
                {"x0": [x1_0, x2_0], "v0": [xdot1_0, xdot2_0]}
            damping_ratio : float | None
                Modal damping ratio ζ (applied equally to all modes)
                for computing damped response.  None → undamped.
        """
        masses   = params.get("masses",  [None, None])
        springs  = params.get("springs", [None, None])
        topology = params.get("topology", "chain")
        ic       = params.get("initial_conditions")
        zeta_in  = params.get("damping_ratio")

        if len(masses) != 2:
            raise ValueError("'masses' must contain exactly 2 entries: [m1, m2].")
        if len(springs) not in (2, 3):
            raise ValueError("'springs' must contain 2 or 3 entries: [k1, k2] or [k1, k2, k3].")

        m1_v, m2_v = masses
        if len(springs) == 2:
            k1_v, k2_v, k3_v = springs[0], springs[1], 0.0
        else:
            k1_v, k2_v, k3_v = springs

        # Topology "free_free" forces k1=k3=0
        if topology == "free_free":
            k1_v, k3_v = 0.0, 0.0

        # ---------------------------------------------------------------
        # Validate completeness for numeric path
        # ---------------------------------------------------------------
        all_numeric = all(
            v is not None for v in [m1_v, m2_v, k1_v, k2_v, k3_v]
        )

        given: dict[str, Any] = {
            "m1": m1_v, "m2": m2_v,
            "k1": k1_v, "k2": k2_v, "k3": k3_v,
            "topology": topology,
        }
        if ic:
            given["x0"] = ic.get("x0", [0, 0])
            given["v0"] = ic.get("v0", [0, 0])
        if zeta_in is not None:
            given["damping_ratio_zeta"] = zeta_in

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # Step 1 – Assemble M and K (symbolic display)
        M_sym, K_sym, m1_s, m2_s, k1_s, k2_s, k3_s = self._build_matrices_sym(
            m1_v, m2_v, k1_v, k2_v, k3_v, topology
        )
        steps.append((
            "Mass matrix M",
            f"M = diag(m1, m2) =\n{sp.sstr(M_sym)}"
        ))
        steps.append((
            "Stiffness matrix K",
            f"K =\n{sp.sstr(K_sym)}"
        ))
        final_answer["M_matrix"] = str(M_sym.tolist())
        final_answer["K_matrix"] = str(K_sym.tolist())

        if not all_numeric:
            steps.append((
                "Note: symbolic parameters detected",
                "Provide numeric values for m1, m2, k1, k2, k3 to obtain "
                "numeric natural frequencies and mode shapes."
            ))
            return SolverResult(
                problem_type="2-DOF Chain – Mass & Stiffness Matrices",
                given=given,
                steps=steps,
                final_answer=final_answer,
                sanity_check="Symbolic mode: only matrices assembled.",
            )

        # ---------------------------------------------------------------
        # Numeric path
        # ---------------------------------------------------------------
        m1, m2 = float(m1_v), float(m2_v)  # type: ignore[arg-type]
        k1, k2, k3 = float(k1_v), float(k2_v), float(k3_v)  # type: ignore[arg-type]

        M_np = np.array([[m1, 0.0], [0.0, m2]])
        K_np = np.array([
            [k1 + k2, -k2],
            [-k2,     k2 + k3],
        ])

        # Step 2 – A = M^{-1/2} K M^{-1/2}
        # For diagonal M: M^{-1/2} = diag(1/sqrt(m1), 1/sqrt(m2))
        M_half_inv = np.diag([1.0 / math.sqrt(m1), 1.0 / math.sqrt(m2)])
        A_np = M_half_inv @ K_np @ M_half_inv

        steps.append((
            "Symmetric eigenvalue matrix A = M^{-1/2} K M^{-1/2}",
            f"A =\n{A_np}"
        ))

        # Step 3 – Eigenvalue problem A u = lambda u
        eigvals, eigvecs = np.linalg.eigh(A_np)   # sorted ascending
        omega_n = np.sqrt(np.maximum(eigvals, 0.0))

        steps.append((
            "Eigenvalues λ (= ω²) and natural frequencies ω_n",
            "\n".join(
                f"  Mode {i+1}: λ_{i+1} = {eigvals[i]:.6f}  →  ω_{i+1} = {omega_n[i]:.6f} rad/s"
                for i in range(2)
            )
        ))

        final_answer["lambda_1"] = round(float(eigvals[0]), 8)
        final_answer["lambda_2"] = round(float(eigvals[1]), 8)
        final_answer["omega_1"]  = round(float(omega_n[0]), 8)
        final_answer["omega_2"]  = round(float(omega_n[1]), 8)

        # Step 4 – Mass-normalised mode shapes U = M^{-1/2} * V_eigvec
        U_cols = M_half_inv @ eigvecs       # columns are mass-normalised vectors
        steps.append((
            "Mass-normalised mode shapes U (columns)",
            f"U =\n{U_cols}\n"
            f"Mode 1: {U_cols[:, 0]}\n"
            f"Mode 2: {U_cols[:, 1]}"
        ))
        final_answer["mode_shape_1"] = U_cols[:, 0].tolist()
        final_answer["mode_shape_2"] = U_cols[:, 1].tolist()

        # Orthogonality check
        UTMU = U_cols.T @ M_np @ U_cols
        UTKU = U_cols.T @ K_np @ U_cols
        steps.append((
            "Orthogonality check: U^T M U and U^T K U",
            f"U^T M U =\n{np.round(UTMU, 10)}\n"
            f"U^T K U =\n{np.round(UTKU, 10)}"
        ))

        # Step 5 – Damping
        if zeta_in is not None:
            zeta = float(zeta_in)
            omega_d = omega_n * math.sqrt(1 - zeta**2)
            sigma   = zeta * omega_n
            steps.append((
                "Damped natural frequencies (modal damping)",
                "\n".join(
                    f"  Mode {i+1}: ζ={zeta:.4f}, "
                    f"σ_{i+1}={sigma[i]:.6f}, "
                    f"ω_d{i+1}={omega_d[i]:.6f} rad/s"
                    for i in range(2)
                )
            ))
            final_answer["zeta"] = zeta
            final_answer["omega_d_1"] = round(float(omega_d[0]), 8)
            final_answer["omega_d_2"] = round(float(omega_d[1]), 8)
            final_answer["sigma_1"]   = round(float(sigma[0]),   8)
            final_answer["sigma_2"]   = round(float(sigma[1]),   8)

        # Step 6 – Initial conditions → modal response
        if ic is not None:
            x0  = np.array(ic.get("x0",  [0.0, 0.0]), dtype=float)
            v0  = np.array(ic.get("v0",  [0.0, 0.0]), dtype=float)

            # Modal initial conditions: eta = U^T M q
            eta0     = U_cols.T @ M_np @ x0
            eta_dot0 = U_cols.T @ M_np @ v0

            steps.append((
                "Modal initial conditions η(0) = U^T M q(0)",
                f"η(0) = {eta0}\n"
                f"η̇(0) = {eta_dot0}"
            ))
            final_answer["modal_ic_eta0"]     = eta0.tolist()
            final_answer["modal_ic_eta_dot0"] = eta_dot0.tolist()

            # Response expression (symbolic — display only)
            if zeta_in is None or float(zeta_in) == 0.0:
                for i in range(2):
                    wn = omega_n[i]
                    A_coeff = eta0[i]
                    B_coeff = eta_dot0[i] / wn if wn > 1e-12 else eta_dot0[i]
                    steps.append((
                        f"Modal response η_{i+1}(t) [undamped]",
                        f"η_{i+1}(t) = {A_coeff:.6f}·cos({wn:.6f}t) "
                        f"+ {B_coeff:.6f}·sin({wn:.6f}t)"
                    ))
            else:
                zeta = float(zeta_in)
                omega_d_arr = omega_n * math.sqrt(1 - zeta**2)
                sigma_arr   = zeta * omega_n
                for i in range(2):
                    wn  = omega_n[i]
                    wd  = omega_d_arr[i]
                    sig = sigma_arr[i]
                    A_c = eta0[i]
                    B_c = (eta_dot0[i] + sig * eta0[i]) / wd if wd > 1e-12 else eta_dot0[i]
                    steps.append((
                        f"Modal response η_{i+1}(t) [damped, ζ={zeta:.4f}]",
                        f"η_{i+1}(t) = e^{{-{sig:.6f}t}}·"
                        f"[{A_c:.6f}·cos({wd:.6f}t) + {B_c:.6f}·sin({wd:.6f}t)]"
                    ))

        # Sanity check
        sanity_parts: list[str] = [
            f"U^T M U diagonal entries: {np.diag(UTMU).round(8).tolist()} (should be [1,1]) "
            + ("✓" if np.allclose(np.diag(UTMU), 1.0, atol=1e-8) else "✗"),
            f"U^T M U off-diagonal: {UTMU[0,1]:.2e} (should be ~0) "
            + ("✓" if abs(UTMU[0, 1]) < 1e-8 else "✗"),
            f"U^T K U diagonal: {np.diag(UTKU).round(8).tolist()} (should be λ values)",
        ]

        return SolverResult(
            problem_type="2-DOF Chain – Modal Analysis",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    def get_input_template(self) -> dict:
        """Return the expected input parameter schema."""
        return {
            "masses": {
                "type": "list[float | None]",
                "description": "[m1, m2] in kg.  None entries → symbolic.",
                "example": [9.0, 9.0],
            },
            "springs": {
                "type": "list[float | None]",
                "description": (
                    "[k1, k2] or [k1, k2, k3] in N/m.\n"
                    "k1 = wall-to-m1 spring,\n"
                    "k2 = m1-to-m2 coupling spring,\n"
                    "k3 = m2-to-wall spring (default 0)."
                ),
                "example": [36.0, 36.0, 36.0],
            },
            "topology": {
                "type": "str",
                "description": (
                    '"chain"              → wall─k1─m1─k2─m2─k3─wall  (default)\n'
                    '"grounded_both_ends" → alias for chain with k1,k3>0\n'
                    '"free_free"          → m1─k2─m2  (k1=k3=0 forced)'
                ),
                "default": "chain",
            },
            "initial_conditions": {
                "type": "dict | None",
                "description": '{"x0": [x1_0, x2_0], "v0": [xdot1_0, xdot2_0]}',
                "example": {"x0": [0.0, 1.0], "v0": [0.0, 0.0]},
            },
            "damping_ratio": {
                "type": "float | None",
                "description": "Modal damping ratio ζ (same for all modes).  None → undamped.",
                "example": 0.15,
            },
        }

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _build_matrices_sym(m1_v, m2_v, k1_v, k2_v, k3_v, topology: str):
        """Build symbolic M and K matrices."""
        m1_s = sp.Symbol("m1") if m1_v is None else sp.nsimplify(m1_v)
        m2_s = sp.Symbol("m2") if m2_v is None else sp.nsimplify(m2_v)
        k1_s = sp.Symbol("k1") if k1_v is None else sp.nsimplify(k1_v)
        k2_s = sp.Symbol("k2") if k2_v is None else sp.nsimplify(k2_v)
        k3_s = sp.Symbol("k3") if k3_v is None else sp.nsimplify(k3_v)

        if topology == "free_free":
            k1_s = sp.Integer(0)
            k3_s = sp.Integer(0)

        M_sym = sp.Matrix([[m1_s, 0], [0, m2_s]])
        K_sym = sp.Matrix([
            [k1_s + k2_s, -k2_s],
            [-k2_s,        k2_s + k3_s],
        ])
        return M_sym, K_sym, m1_s, m2_s, k1_s, k2_s, k3_s
