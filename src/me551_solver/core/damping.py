"""Proportional Damping 솔버.

비례감쇠(C = αM + βK) 시스템의 모달 분해 및 감쇠 응답을 계산한다.
비비례감쇠의 경우 state-space 접근을 안내한다.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.linalg import eigh

from .base import BaseSolver, SolverResult
from .modal import ModalSolver


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _format_matrix(mat: np.ndarray, name: str = "") -> str:
    """Pretty-format a numpy matrix for display."""
    lines = []
    if name:
        lines.append(f"{name} =")
    for row in mat:
        row_str = "  [" + ", ".join(f"{x: .6g}" for x in row) + "]"
        lines.append(row_str)
    return "\n".join(lines)


def _format_vector(vec: np.ndarray, name: str = "") -> str:
    """Pretty-format a numpy vector for display."""
    vec_str = "[" + ", ".join(f"{x: .6g}" for x in vec) + "]"
    if name:
        return f"{name} = {vec_str}"
    return vec_str


def _check_proportional_damping(M: np.ndarray, K: np.ndarray,
                                 C: np.ndarray) -> dict[str, Any]:
    """Check if C = αM + βK for some α, β.

    Method 1: Solve the overdetermined system C = αM + βK using least-squares.
    Method 2: Check commutativity M⁻¹C M⁻¹K = M⁻¹K M⁻¹C.

    Returns:
        dict with keys:
            is_proportional: bool
            alpha: float or None
            beta: float or None
            residual: float (Frobenius norm of C - αM - βK)
            commutativity_check: bool
    """
    n = M.shape[0]

    # Flatten matrices for least-squares: vec(C) = α vec(M) + β vec(K)
    M_flat = M.flatten()
    K_flat = K.flatten()
    C_flat = C.flatten()

    # Build the system: A @ [alpha, beta]^T = C_flat
    A_mat = np.column_stack([M_flat, K_flat])
    result_ls = np.linalg.lstsq(A_mat, C_flat, rcond=None)
    coeffs = result_ls[0]
    alpha, beta = coeffs[0], coeffs[1]

    # Compute residual
    C_reconstructed = alpha * M + beta * K
    residual = np.linalg.norm(C - C_reconstructed, 'fro')
    C_norm = np.linalg.norm(C, 'fro')
    relative_residual = residual / C_norm if C_norm > 1e-15 else residual

    # Commutativity check: M⁻¹C M⁻¹K == M⁻¹K M⁻¹C
    M_inv = np.linalg.inv(M)
    MiC = M_inv @ C
    MiK = M_inv @ K
    comm_diff = MiC @ MiK - MiK @ MiC
    comm_norm = np.linalg.norm(comm_diff, 'fro')
    is_commutative = comm_norm < 1e-8

    # Proportionality threshold
    is_proportional = relative_residual < 1e-8 and is_commutative

    return {
        "is_proportional": is_proportional,
        "alpha": float(alpha),
        "beta": float(beta),
        "residual": float(residual),
        "relative_residual": float(relative_residual),
        "commutativity_check": is_commutative,
        "commutativity_norm": float(comm_norm),
    }


def _solve_nonproportional_state_space(M: np.ndarray, K: np.ndarray,
                                        C: np.ndarray) -> dict[str, Any]:
    """Solve nonproportional damping via state-space eigenvalue analysis.

    State-space form:
        [M  0] {q̈}   [C  K] {q̇}
        [0 -K] {q̇} + [K  0] {q}  = 0

    Or equivalently, first-order form:
        ẋ = A x,   x = [q, q̇]^T
        A = [[0, I], [-M⁻¹K, -M⁻¹C]]

    Returns complex eigenvalues from which ω and ζ are extracted.
    """
    n = M.shape[0]
    M_inv = np.linalg.inv(M)

    # Build state matrix A (2n x 2n)
    A = np.zeros((2 * n, 2 * n))
    A[:n, n:] = np.eye(n)
    A[n:, :n] = -M_inv @ K
    A[n:, n:] = -M_inv @ C

    # Compute eigenvalues
    eigvals = np.linalg.eigvals(A)

    # Extract complex conjugate pairs (keep those with positive imaginary part)
    # Sort by imaginary part magnitude
    pairs = []
    used = set()
    for i, ev in enumerate(eigvals):
        if i in used:
            continue
        if abs(ev.imag) > 1e-10:
            # Find conjugate pair
            for j, ev2 in enumerate(eigvals):
                if j > i and j not in used:
                    if abs(ev + ev2.conjugate()) < 1e-8 or abs(ev.real - ev2.real) < 1e-8:
                        if ev.imag > 0:
                            pairs.append(ev)
                        else:
                            pairs.append(ev2)
                        used.add(i)
                        used.add(j)
                        break
        elif abs(ev.imag) <= 1e-10:
            # Real eigenvalue (overdamped)
            pairs.append(ev)
            used.add(i)

    # Sort by natural frequency magnitude
    pairs.sort(key=lambda x: abs(x))

    # Remove duplicate real eigenvalues (they come in pairs for overdamped)
    # Keep unique pairs
    natural_frequencies = []
    damping_ratios = []
    damped_frequencies = []

    processed = set()
    for ev in pairs:
        sigma = -ev.real   # decay rate
        omega_d = abs(ev.imag)   # damped frequency
        omega_n = abs(ev)        # natural frequency = |λ|

        if omega_n < 1e-14:
            continue

        # Avoid counting the same mode twice
        key = round(omega_n, 6)
        if key in processed:
            continue
        processed.add(key)

        zeta = sigma / omega_n

        natural_frequencies.append(omega_n)
        damping_ratios.append(zeta)
        damped_frequencies.append(omega_d)

    return {
        "natural_frequencies": natural_frequencies,
        "damping_ratios": damping_ratios,
        "damped_frequencies": damped_frequencies,
        "state_matrix_eigenvalues": eigvals.tolist(),
    }


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class DampingSolver(BaseSolver):
    """비례감쇠(Proportional Damping) 및 비비례감쇠 시스템의 응답 해석.

    Proportional damping: C = αM + βK → decoupled modal equations.
    Nonproportional damping: state-space eigenvalue analysis.
    """

    def solve(self, params: dict) -> SolverResult:
        """Perform damped modal analysis.

        Parameters:
            params: dict with keys:
                M: mass matrix (list of lists or numpy array)
                K: stiffness matrix (list of lists or numpy array)
                C: damping matrix (optional, list of lists or numpy array)
                damping_type: "proportional", "nonproportional", or "check"
                alpha: Rayleigh damping coefficient α (optional)
                beta: Rayleigh damping coefficient β (optional)
                zeta: list of modal damping ratios (optional, alternative to C)
                initial_q: initial displacement vector (optional)
                initial_qdot: initial velocity vector (optional)
        """
        # ---------------------------------------------------------------
        # 0. Parse inputs
        # ---------------------------------------------------------------
        M_input = params.get("M")
        K_input = params.get("K")
        if M_input is None or K_input is None:
            raise ValueError("Mass matrix 'M' and stiffness matrix 'K' are required.")

        M = np.array(M_input, dtype=float)
        K = np.array(K_input, dtype=float)
        n_dof = M.shape[0]

        C_input = params.get("C")
        C = np.array(C_input, dtype=float) if C_input is not None else None

        alpha = params.get("alpha")
        beta = params.get("beta")
        zeta_list = params.get("zeta")
        damping_type = params.get("damping_type", "check")

        initial_q = params.get("initial_q")
        initial_qdot = params.get("initial_qdot")
        q0 = np.array(initial_q, dtype=float) if initial_q is not None else np.zeros(n_dof)
        qdot0 = np.array(initial_qdot, dtype=float) if initial_qdot is not None else np.zeros(n_dof)

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # ---------------------------------------------------------------
        # 1. Perform undamped modal analysis first
        # ---------------------------------------------------------------
        modal_solver = ModalSolver()
        modal_result = modal_solver.solve({
            "M": M.tolist(),
            "K": K.tolist(),
            "initial_q": q0.tolist(),
            "initial_qdot": qdot0.tolist(),
        })

        eigenvalues = np.array(modal_result.final_answer["eigenvalues"])
        omega = np.array(modal_result.final_answer["natural_frequencies"])
        U = np.array(modal_result.final_answer["modal_matrix"])

        steps.append((
            "Undamped modal analysis (via ModalSolver)",
            f"Eigenvalues: {eigenvalues}\n"
            f"Natural frequencies: {omega} rad/s\n"
            f"Modal matrix U:\n{_format_matrix(U, 'U')}"
        ))
        final_answer["eigenvalues"] = eigenvalues.tolist()
        final_answer["undamped_natural_frequencies"] = omega.tolist()

        # ---------------------------------------------------------------
        # 2. Build or check damping matrix C
        # ---------------------------------------------------------------
        if C is None and alpha is not None and beta is not None:
            # Build C from Rayleigh coefficients
            C = alpha * M + beta * K
            steps.append((
                "Construct damping matrix: C = αM + βK",
                f"α = {alpha}, β = {beta}\n"
                + _format_matrix(C, "C")
            ))
            damping_type = "proportional"
        elif C is None and zeta_list is not None:
            # Damping ratios given directly
            zeta_arr = np.array(zeta_list, dtype=float)
            steps.append((
                "Modal damping ratios given directly",
                _format_vector(zeta_arr, "ζ")
            ))
            damping_type = "proportional_modal"
        elif C is not None:
            steps.append((
                "Damping matrix C",
                _format_matrix(C, "C")
            ))
        else:
            raise ValueError(
                "Either 'C' (damping matrix), 'alpha'+'beta' (Rayleigh coefficients), "
                "or 'zeta' (modal damping ratios) must be provided."
            )

        # ---------------------------------------------------------------
        # 3. Check proportional damping
        # ---------------------------------------------------------------
        is_proportional = False
        alpha_val = alpha
        beta_val = beta

        if damping_type == "proportional_modal":
            is_proportional = True
            zeta_arr = np.array(zeta_list, dtype=float)
            steps.append((
                "Damping type: proportional (modal damping ratios specified)",
                "Modal damping ratios directly applied to each mode."
            ))
        elif damping_type == "proportional":
            is_proportional = True
            if alpha_val is not None and beta_val is not None:
                steps.append((
                    "Check proportional damping: C = αM + βK",
                    f"Proportional damping confirmed with α = {alpha_val}, β = {beta_val}"
                ))
            else:
                # Verify proportionality
                prop_check = _check_proportional_damping(M, K, C)
                is_proportional = prop_check["is_proportional"]
                alpha_val = prop_check["alpha"]
                beta_val = prop_check["beta"]
                steps.append((
                    "Check proportional damping: C = αM + βK",
                    f"Least-squares fit: α = {alpha_val:.6g}, β = {beta_val:.6g}\n"
                    f"Residual (Frobenius): {prop_check['residual']:.2e}\n"
                    f"Commutativity M⁻¹C M⁻¹K = M⁻¹K M⁻¹C: "
                    f"{'PASS ✓' if prop_check['commutativity_check'] else 'FAIL ✗'}\n"
                    f"Result: {'PROPORTIONAL ✓' if is_proportional else 'NOT PROPORTIONAL ✗'}"
                ))
        elif damping_type == "check" and C is not None:
            prop_check = _check_proportional_damping(M, K, C)
            is_proportional = prop_check["is_proportional"]
            alpha_val = prop_check["alpha"]
            beta_val = prop_check["beta"]
            steps.append((
                "Check proportional damping: C = αM + βK",
                f"Least-squares fit: α = {alpha_val:.6g}, β = {beta_val:.6g}\n"
                f"Residual (Frobenius): {prop_check['residual']:.2e}\n"
                f"Relative residual: {prop_check['relative_residual']:.2e}\n"
                f"Commutativity M⁻¹C M⁻¹K = M⁻¹K M⁻¹C: "
                f"{'PASS ✓' if prop_check['commutativity_check'] else 'FAIL ✗'}\n"
                f"Commutativity norm: {prop_check['commutativity_norm']:.2e}\n"
                f"Result: {'PROPORTIONAL ✓' if is_proportional else 'NOT PROPORTIONAL ✗'}"
            ))
            final_answer["is_proportional"] = is_proportional
        elif damping_type == "nonproportional":
            is_proportional = False
            steps.append((
                "Damping type: nonproportional (specified by user)",
                "Skipping proportionality check. Will use state-space approach."
            ))

        final_answer["is_proportional"] = is_proportional
        if alpha_val is not None:
            final_answer["alpha"] = float(alpha_val)
        if beta_val is not None:
            final_answer["beta"] = float(beta_val)

        # ---------------------------------------------------------------
        # 4. Compute modal damping ratios
        # ---------------------------------------------------------------
        if is_proportional:
            if damping_type == "proportional_modal":
                zeta = np.array(zeta_list, dtype=float)
            else:
                # ζ_r = (α + β ω_r²) / (2 ω_r)
                zeta = np.zeros(n_dof)
                for r in range(n_dof):
                    if omega[r] > 1e-14:
                        zeta[r] = (alpha_val + beta_val * omega[r] ** 2) / (2 * omega[r])
                    else:
                        zeta[r] = 0.0

            zeta_strs = []
            for r in range(n_dof):
                if damping_type != "proportional_modal" and omega[r] > 1e-14:
                    zeta_strs.append(
                        f"2ζ_{r+1} ω_{r+1} = α + β ω_{r+1}² = "
                        f"{alpha_val:.6g} + {beta_val:.6g} × {omega[r]**2:.6g} = "
                        f"{alpha_val + beta_val * omega[r]**2:.6g}\n"
                        f"  → ζ_{r+1} = {zeta[r]:.6g}"
                    )
                else:
                    zeta_strs.append(f"ζ_{r+1} = {zeta[r]:.6g}")

            steps.append((
                "Modal damping: 2ζ_r ω_r = α + β ω_r²",
                "\n".join(zeta_strs)
            ))

            for r in range(n_dof):
                final_answer[f"zeta_{r+1}"] = float(zeta[r])

            # ---------------------------------------------------------------
            # 5. Damped natural frequencies
            # ---------------------------------------------------------------
            omega_d = np.zeros(n_dof)
            omega_d_strs = []
            for r in range(n_dof):
                if zeta[r] < 1.0 and omega[r] > 1e-14:
                    omega_d[r] = omega[r] * np.sqrt(1 - zeta[r] ** 2)
                    omega_d_strs.append(
                        f"ω_d,{r+1} = ω_{r+1} √(1 - ζ_{r+1}²) = "
                        f"{omega[r]:.6g} × √(1 - {zeta[r]:.6g}²) = {omega_d[r]:.6g} rad/s"
                    )
                elif zeta[r] >= 1.0:
                    omega_d[r] = 0.0
                    omega_d_strs.append(
                        f"ω_d,{r+1} = 0 (overdamped, ζ_{r+1} = {zeta[r]:.6g} ≥ 1)"
                    )
                else:
                    omega_d[r] = 0.0
                    omega_d_strs.append(
                        f"ω_d,{r+1} = 0 (rigid-body mode)"
                    )

            steps.append((
                "Damped natural frequencies: ω_d,r = ω_r √(1 - ζ_r²)",
                "\n".join(omega_d_strs)
            ))

            for r in range(n_dof):
                final_answer[f"omega_{r+1}"] = float(omega[r])
                final_answer[f"omega_d_{r+1}"] = float(omega_d[r])

            # ---------------------------------------------------------------
            # 6. Modal equations
            # ---------------------------------------------------------------
            modal_eq_strs = []
            for r in range(n_dof):
                modal_eq_strs.append(
                    f"η̈_{r+1} + 2×{zeta[r]:.6g}×{omega[r]:.6g} η̇_{r+1} "
                    f"+ {omega[r]**2:.6g} η_{r+1} = 0"
                )
            steps.append((
                "Modal equations: η̈_r + 2ζ_r ω_r η̇_r + ω_r² η_r = 0",
                "\n".join(modal_eq_strs)
            ))

            # ---------------------------------------------------------------
            # 7. Modal response (damped free vibration)
            # ---------------------------------------------------------------
            eta0 = np.array(modal_result.final_answer["eta_0"])
            etadot0 = np.array(modal_result.final_answer["etadot_0"])

            modal_resp_strs = []
            A_coeffs = np.zeros(n_dof)
            B_coeffs = np.zeros(n_dof)

            for r in range(n_dof):
                sigma_r = zeta[r] * omega[r]
                if omega_d[r] > 1e-14:
                    # Underdamped case
                    A_r = eta0[r]
                    B_r = (etadot0[r] + sigma_r * eta0[r]) / omega_d[r]
                    A_coeffs[r] = A_r
                    B_coeffs[r] = B_r
                    modal_resp_strs.append(
                        f"η_{r+1}(t) = e^(-{sigma_r:.6g} t) "
                        f"[{A_r:.6g} cos({omega_d[r]:.6g} t) "
                        f"+ {B_r:.6g} sin({omega_d[r]:.6g} t)]"
                    )
                elif omega[r] > 1e-14 and zeta[r] >= 1.0:
                    # Overdamped (critical or overdamped)
                    modal_resp_strs.append(
                        f"η_{r+1}(t) = overdamped response (ζ_{r+1} ≥ 1)"
                    )
                else:
                    # Rigid-body / zero frequency
                    A_coeffs[r] = eta0[r]
                    B_coeffs[r] = etadot0[r]
                    modal_resp_strs.append(
                        f"η_{r+1}(t) = {eta0[r]:.6g} + {etadot0[r]:.6g} t"
                    )

            steps.append((
                "Modal response: η_r(t) = e^(-ζ_r ω_r t)[A_r cos(ω_d,r t) + B_r sin(ω_d,r t)]",
                "\n".join(modal_resp_strs)
            ))

            # ---------------------------------------------------------------
            # 8. Physical response
            # ---------------------------------------------------------------
            phys_resp_strs = []
            for i in range(n_dof):
                terms = []
                for r in range(n_dof):
                    u_ir = U[i, r]
                    if abs(u_ir) < 1e-14:
                        continue
                    sigma_r = zeta[r] * omega[r]
                    if omega_d[r] > 1e-14:
                        a_phys = u_ir * A_coeffs[r]
                        b_phys = u_ir * B_coeffs[r]
                        terms.append(
                            f"e^(-{sigma_r:.4g}t)"
                            f"[{a_phys:.6g} cos({omega_d[r]:.4g}t)"
                            f" + {b_phys:.6g} sin({omega_d[r]:.4g}t)]"
                        )
                    else:
                        terms.append(f"{u_ir * A_coeffs[r]:.6g}")

                if terms:
                    phys_resp_strs.append(f"q_{i+1}(t) = " + " + ".join(terms))
                else:
                    phys_resp_strs.append(f"q_{i+1}(t) = 0")

            steps.append((
                "Physical response: q(t) = U η(t)",
                "\n".join(phys_resp_strs)
            ))

        # ---------------------------------------------------------------
        # Nonproportional damping: state-space approach
        # ---------------------------------------------------------------
        else:
            if C is not None:
                steps.append((
                    "Nonproportional damping detected",
                    "C ≠ αM + βK for any α, β.\n"
                    "Using state-space eigenvalue analysis to extract\n"
                    "complex natural frequencies and damping ratios."
                ))

                ss_result = _solve_nonproportional_state_space(M, K, C)

                omega_np = ss_result["natural_frequencies"]
                zeta_np = ss_result["damping_ratios"]
                omega_d_np = ss_result["damped_frequencies"]

                ss_strs = []
                for r in range(len(omega_np)):
                    ss_strs.append(
                        f"Mode {r+1}: ω_{r+1} = {omega_np[r]:.6g} rad/s, "
                        f"ζ_{r+1} = {zeta_np[r]:.6g}, "
                        f"ω_d,{r+1} = {omega_d_np[r]:.6g} rad/s"
                    )

                steps.append((
                    "State-space eigenvalue results",
                    "State matrix A = [[0, I], [-M⁻¹K, -M⁻¹C]]\n"
                    "Complex eigenvalues: λ = -σ ± jω_d\n"
                    "  → ω_n = |λ|, ζ = σ/ω_n, ω_d = ω_n√(1-ζ²)\n\n"
                    + "\n".join(ss_strs)
                ))

                for r in range(len(omega_np)):
                    final_answer[f"omega_{r+1}"] = omega_np[r]
                    final_answer[f"zeta_{r+1}"] = zeta_np[r]
                    final_answer[f"omega_d_{r+1}"] = omega_d_np[r]

                steps.append((
                    "Note on nonproportional damping",
                    "For nonproportional damping, the undamped mode shapes do NOT\n"
                    "decouple the equations of motion. The complex mode shapes from\n"
                    "the state-space analysis are required for the full response.\n"
                    "A full complex-mode response is not yet implemented in this version."
                ))
            else:
                steps.append((
                    "Nonproportional damping",
                    "No damping matrix C provided for nonproportional analysis.\n"
                    "Please provide the full damping matrix C."
                ))

        # ---------------------------------------------------------------
        # Sanity check
        # ---------------------------------------------------------------
        sanity_parts: list[str] = []

        if is_proportional:
            sanity_parts.append("Damping type: PROPORTIONAL (C = αM + βK)")
            # Check all damping ratios are non-negative
            all_positive = all(zeta[r] >= 0 for r in range(n_dof))
            sanity_parts.append(
                f"All damping ratios non-negative: "
                f"{'YES ✓' if all_positive else 'WARNING: negative damping detected'}"
            )
            # Check underdamped
            all_underdamped = all(zeta[r] < 1.0 for r in range(n_dof))
            sanity_parts.append(
                f"All modes underdamped (ζ < 1): "
                f"{'YES ✓' if all_underdamped else 'NO — some modes are critically/overdamped'}"
            )
        else:
            sanity_parts.append("Damping type: NONPROPORTIONAL")
            sanity_parts.append(
                "State-space approach used for eigenvalue extraction."
            )

        sanity_check_str = "\n".join(sanity_parts)

        return SolverResult(
            problem_type="Damped MDOF Modal Analysis",
            given={
                "M": M.tolist(),
                "K": K.tolist(),
                "C": C.tolist() if C is not None else None,
                "damping_type": damping_type,
                "alpha": alpha,
                "beta": beta,
                "initial_q": q0.tolist(),
                "initial_qdot": qdot0.tolist(),
            },
            steps=steps,
            final_answer=final_answer,
            sanity_check=sanity_check_str,
        )

    def get_input_template(self) -> dict:
        """필요한 입력 파라미터 설명을 반환한다."""
        return {
            "M": {
                "type": "list[list[float]] or np.ndarray",
                "required": True,
                "description": "Mass matrix (n x n, symmetric positive-definite).",
                "example": [[9, 0], [0, 9]],
            },
            "K": {
                "type": "list[list[float]] or np.ndarray",
                "required": True,
                "description": "Stiffness matrix (n x n, symmetric positive semi-definite).",
                "example": [[72, -36], [-36, 72]],
            },
            "C": {
                "type": "list[list[float]] or np.ndarray",
                "required": False,
                "description": (
                    "Damping matrix (n x n). Required unless alpha/beta or zeta provided."
                ),
                "example": [[2, -1], [-1, 2]],
            },
            "damping_type": {
                "type": "str",
                "required": False,
                "default": "check",
                "description": (
                    "'proportional', 'nonproportional', or 'check' (auto-detect)."
                ),
            },
            "alpha": {
                "type": "float",
                "required": False,
                "description": "Rayleigh damping mass coefficient α (C = αM + βK).",
            },
            "beta": {
                "type": "float",
                "required": False,
                "description": "Rayleigh damping stiffness coefficient β (C = αM + βK).",
            },
            "zeta": {
                "type": "list[float]",
                "required": False,
                "description": (
                    "Modal damping ratios ζ_r for each mode. "
                    "Alternative to specifying C or α, β."
                ),
                "example": [0.05, 0.1],
            },
            "initial_q": {
                "type": "list[float]",
                "required": False,
                "default": "zeros",
                "description": "Initial displacement vector q(0).",
            },
            "initial_qdot": {
                "type": "list[float]",
                "required": False,
                "default": "zeros",
                "description": "Initial velocity vector q̇(0).",
            },
        }
