"""N-DOF Modal Analysis 솔버 (undamped).

Generalized eigenvalue problem K u = λ M u 을 풀고,
질량 정규화 모드형상, 직교성 검증, 모달/물리 응답을 계산한다.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.linalg import eigh

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _fix_sign_convention(vectors: np.ndarray) -> np.ndarray:
    """Fix eigenvector sign convention: first nonzero element positive.

    Parameters:
        vectors: 2-D array where each column is an eigenvector.

    Returns:
        Copy with consistent sign convention applied per column.
    """
    result = vectors.copy()
    for j in range(result.shape[1]):
        col = result[:, j]
        for val in col:
            if abs(val) > 1e-12:
                if val < 0:
                    result[:, j] = -col
                break
    return result


def _mass_normalize(vectors: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Mass-normalize eigenvectors so that u_i^T M u_i = 1.

    Parameters:
        vectors: 2-D array where each column is an eigenvector.
        M: mass matrix (n x n).

    Returns:
        Mass-normalized eigenvector matrix (each column normalized).
    """
    result = vectors.copy()
    n_modes = result.shape[1]
    for j in range(n_modes):
        v = result[:, j]
        scale = np.sqrt(v @ M @ v)
        result[:, j] = v / scale
    return result


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


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class ModalSolver(BaseSolver):
    """N-DOF 시스템의 모달 분석 – 고유진동수, 모드형상, 모달 좌표 응답.

    Uses scipy.linalg.eigh for the generalized symmetric eigenvalue problem
    K u = λ M u, which is the most reliable numerical method.
    """

    def solve(self, params: dict) -> SolverResult:
        """Perform undamped modal analysis.

        Parameters:
            params: dict with keys:
                M: mass matrix (list of lists or numpy array)
                K: stiffness matrix (list of lists or numpy array)
                initial_q: initial displacement vector (optional)
                initial_qdot: initial velocity vector (optional)
                compute: list of computation stages to include (optional)
                    defaults to all stages
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

        initial_q = params.get("initial_q")
        initial_qdot = params.get("initial_qdot")
        if initial_q is not None:
            q0 = np.array(initial_q, dtype=float)
        else:
            q0 = np.zeros(n_dof)
        if initial_qdot is not None:
            qdot0 = np.array(initial_qdot, dtype=float)
        else:
            qdot0 = np.zeros(n_dof)

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # ---------------------------------------------------------------
        # 1. Display M and K
        # ---------------------------------------------------------------
        steps.append((
            "Mass matrix M",
            _format_matrix(M, "M")
        ))
        steps.append((
            "Stiffness matrix K",
            _format_matrix(K, "K")
        ))

        # ---------------------------------------------------------------
        # 2. Setup eigenvalue problem
        # ---------------------------------------------------------------
        steps.append((
            "Generalized eigenvalue problem: K u = λ M u",
            "Solve for eigenvalues λ and eigenvectors u using the\n"
            "generalized symmetric eigenvalue problem."
        ))

        # ---------------------------------------------------------------
        # 3. Solve using scipy.linalg.eigh
        # ---------------------------------------------------------------
        steps.append((
            "Convert to standard form: M⁻¹K v = λ v (or use scipy.linalg.eigh)",
            "Using scipy.linalg.eigh(K, M) for numerically stable solution\n"
            "of the generalized symmetric eigenvalue problem."
        ))

        eigenvalues, raw_eigvecs = eigh(K, M)
        # eigh returns eigenvalues in ascending order and corresponding eigenvectors

        # ---------------------------------------------------------------
        # 4. Eigenvalues
        # ---------------------------------------------------------------
        # Sort (eigh already sorts, but be explicit)
        sort_idx = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[sort_idx]
        raw_eigvecs = raw_eigvecs[:, sort_idx]

        eigenval_strs = []
        for i, lam in enumerate(eigenvalues):
            eigenval_strs.append(f"λ_{i+1} = {lam:.6g}")
        steps.append((
            "Eigenvalues λ_i",
            "\n".join(eigenval_strs)
        ))
        final_answer["eigenvalues"] = eigenvalues.tolist()

        # ---------------------------------------------------------------
        # 5. Natural frequencies
        # ---------------------------------------------------------------
        # Handle potential tiny negative eigenvalues from numerical noise
        omega = np.sqrt(np.maximum(eigenvalues, 0.0))
        omega_strs = []
        for i, w in enumerate(omega):
            omega_strs.append(f"ω_{i+1} = √{eigenvalues[i]:.6g} = {w:.6g} rad/s")
        steps.append((
            "Natural frequencies ω_i = √λ_i",
            "\n".join(omega_strs)
        ))
        final_answer["natural_frequencies"] = omega.tolist()

        # ---------------------------------------------------------------
        # 6. Raw eigenvectors
        # ---------------------------------------------------------------
        raw_eigvecs = _fix_sign_convention(raw_eigvecs)
        eigvec_strs = []
        for i in range(n_dof):
            eigvec_strs.append(
                f"v_{i+1} = " + _format_vector(raw_eigvecs[:, i])
            )
        steps.append((
            "Raw eigenvectors",
            "\n".join(eigvec_strs)
        ))

        # ---------------------------------------------------------------
        # 7. Mass-normalized eigenvectors
        # ---------------------------------------------------------------
        U = _mass_normalize(raw_eigvecs, M)
        U = _fix_sign_convention(U)
        norm_strs = []
        for i in range(n_dof):
            check_val = U[:, i] @ M @ U[:, i]
            norm_strs.append(
                f"u_{i+1} = {_format_vector(U[:, i])}"
                f"  (u_{i+1}ᵀ M u_{i+1} = {check_val:.6g})"
            )
        steps.append((
            "Mass-normalized eigenvectors u_i (u_iᵀ M u_i = 1)",
            "\n".join(norm_strs)
        ))
        final_answer["mass_normalized_eigenvectors"] = U.tolist()

        # ---------------------------------------------------------------
        # 8. Modal matrix U
        # ---------------------------------------------------------------
        steps.append((
            "Modal matrix U = [u_1, u_2, ...]",
            _format_matrix(U, "U")
        ))
        final_answer["modal_matrix"] = U.tolist()

        # ---------------------------------------------------------------
        # 9. Orthogonality check: U^T M U = I
        # ---------------------------------------------------------------
        UtMU = U.T @ M @ U
        identity_check = np.allclose(UtMU, np.eye(n_dof), atol=1e-10)
        steps.append((
            "Orthogonality check: Uᵀ M U = I",
            _format_matrix(UtMU, "Uᵀ M U") +
            f"\n\nIdentity check: {'PASS ✓' if identity_check else 'FAIL ✗'}"
        ))

        # ---------------------------------------------------------------
        # 10. Orthogonality check: U^T K U = Λ
        # ---------------------------------------------------------------
        UtKU = U.T @ K @ U
        Lambda_diag = np.diag(eigenvalues)
        spectral_check = np.allclose(UtKU, Lambda_diag, atol=1e-8)
        steps.append((
            "Orthogonality check: Uᵀ K U = Λ",
            _format_matrix(UtKU, "Uᵀ K U") +
            f"\n\nExpected diagonal: {eigenvalues}"
            f"\nSpectral check: {'PASS ✓' if spectral_check else 'FAIL ✗'}"
        ))

        # ---------------------------------------------------------------
        # 11. Modal initial conditions
        # ---------------------------------------------------------------
        eta0 = U.T @ M @ q0
        etadot0 = U.T @ M @ qdot0
        ic_strs = []
        for i in range(n_dof):
            ic_strs.append(
                f"η_{i+1}(0) = u_{i+1}ᵀ M q(0) = {eta0[i]:.6g}"
            )
            ic_strs.append(
                f"η̇_{i+1}(0) = u_{i+1}ᵀ M q̇(0) = {etadot0[i]:.6g}"
            )
        steps.append((
            "Modal initial conditions: η(0) = Uᵀ M q(0), η̇(0) = Uᵀ M q̇(0)",
            "\n".join(ic_strs)
        ))
        final_answer["eta_0"] = eta0.tolist()
        final_answer["etadot_0"] = etadot0.tolist()

        # ---------------------------------------------------------------
        # 12. Modal response expressions
        # ---------------------------------------------------------------
        modal_resp_strs = []
        for i in range(n_dof):
            if abs(omega[i]) < 1e-14:
                # Rigid-body mode
                modal_resp_strs.append(
                    f"η_{i+1}(t) = {eta0[i]:.6g} + {etadot0[i]:.6g} t  (rigid-body mode)"
                )
            else:
                A_coeff = eta0[i]
                B_coeff = etadot0[i] / omega[i]
                modal_resp_strs.append(
                    f"η_{i+1}(t) = {A_coeff:.6g} cos({omega[i]:.6g} t)"
                    f" + {B_coeff:.6g} sin({omega[i]:.6g} t)"
                )
        steps.append((
            "Modal response: η_r(t) = η_r(0)cos(ω_r t) + η̇_r(0)/ω_r sin(ω_r t)",
            "\n".join(modal_resp_strs)
        ))

        # ---------------------------------------------------------------
        # 13. Physical response expression
        # ---------------------------------------------------------------
        steps.append((
            "Physical response: q(t) = U η(t) = Σ u_r η_r(t)",
            "Each physical coordinate q_i(t) is obtained by summing\n"
            "the contributions from all modal coordinates:\n"
            "  q_i(t) = Σ_r U[i,r] · η_r(t)"
        ))

        # ---------------------------------------------------------------
        # 14. Final physical response expressions
        # ---------------------------------------------------------------
        phys_resp_strs = []
        for i in range(n_dof):
            terms = []
            for r in range(n_dof):
                u_ir = U[i, r]
                if abs(u_ir) < 1e-14:
                    continue

                if abs(omega[r]) < 1e-14:
                    # Rigid-body mode
                    A_coeff = eta0[r]
                    B_coeff = etadot0[r]
                    if abs(A_coeff * u_ir) > 1e-14:
                        terms.append(f"{u_ir * A_coeff:.6g}")
                    if abs(B_coeff * u_ir) > 1e-14:
                        terms.append(f"{u_ir * B_coeff:.6g} t")
                else:
                    A_coeff = u_ir * eta0[r]
                    B_coeff = u_ir * etadot0[r] / omega[r]
                    cos_term = f"{A_coeff:.6g} cos({omega[r]:.6g} t)" if abs(A_coeff) > 1e-14 else ""
                    sin_term = f"{B_coeff:.6g} sin({omega[r]:.6g} t)" if abs(B_coeff) > 1e-14 else ""
                    if cos_term and sin_term:
                        if B_coeff >= 0:
                            terms.append(f"{cos_term} + {sin_term}")
                        else:
                            terms.append(f"{cos_term} - {abs(B_coeff):.6g} sin({omega[r]:.6g} t)")
                    elif cos_term:
                        terms.append(cos_term)
                    elif sin_term:
                        terms.append(sin_term)

            if terms:
                phys_resp_strs.append(f"q_{i+1}(t) = " + " + ".join(terms))
            else:
                phys_resp_strs.append(f"q_{i+1}(t) = 0")

        steps.append((
            "Final physical response expressions",
            "\n".join(phys_resp_strs)
        ))

        # Store omega values with numbered keys for easy access
        for i in range(n_dof):
            final_answer[f"omega_{i+1}"] = float(omega[i])
            final_answer[f"lambda_{i+1}"] = float(eigenvalues[i])

        # ---------------------------------------------------------------
        # Evaluation at specific time (optional)
        # ---------------------------------------------------------------
        t_eval = params.get("t_eval")
        if t_eval is not None:
            t_arr = np.atleast_1d(np.array(t_eval, dtype=float))
            for t_val in t_arr:
                # Compute modal response at time t
                eta_t = np.zeros(n_dof)
                for r in range(n_dof):
                    if abs(omega[r]) < 1e-14:
                        eta_t[r] = eta0[r] + etadot0[r] * t_val
                    else:
                        eta_t[r] = (eta0[r] * np.cos(omega[r] * t_val)
                                    + etadot0[r] / omega[r] * np.sin(omega[r] * t_val))
                q_t = U @ eta_t
                for i in range(n_dof):
                    final_answer[f"q{i+1}_at_t{t_val}"] = float(q_t[i])

        # ---------------------------------------------------------------
        # Sanity check
        # ---------------------------------------------------------------
        sanity_parts: list[str] = []
        # All eigenvalues should be non-negative for a stable undamped system
        if np.all(eigenvalues >= -1e-10):
            sanity_parts.append(
                "All eigenvalues are non-negative → system is stable (positive semi-definite)"
            )
        else:
            sanity_parts.append(
                "WARNING: Negative eigenvalue(s) detected → system may be unstable"
            )

        # Mass-orthonormality
        if identity_check:
            sanity_parts.append("Uᵀ M U = I verified ✓")
        else:
            sanity_parts.append("WARNING: Uᵀ M U ≠ I — mass orthonormality failed")

        # Spectral decomposition
        if spectral_check:
            sanity_parts.append("Uᵀ K U = Λ verified ✓")
        else:
            sanity_parts.append("WARNING: Uᵀ K U ≠ Λ — spectral decomposition check failed")

        # Initial condition reconstruction
        q0_reconstructed = U @ eta0
        if np.allclose(q0_reconstructed, q0, atol=1e-8):
            sanity_parts.append(
                f"IC reconstruction: U η(0) = q(0) verified ✓"
            )
        else:
            sanity_parts.append(
                f"WARNING: IC reconstruction mismatch: "
                f"U η(0) = {q0_reconstructed} ≠ q(0) = {q0}"
            )

        sanity_check_str = "\n".join(sanity_parts)

        return SolverResult(
            problem_type="Modal Analysis (Undamped MDOF)",
            given={
                "M": M.tolist(),
                "K": K.tolist(),
                "n_dof": n_dof,
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
            "initial_q": {
                "type": "list[float]",
                "required": False,
                "default": "zeros",
                "description": "Initial displacement vector q(0).",
                "example": [1, 1],
            },
            "initial_qdot": {
                "type": "list[float]",
                "required": False,
                "default": "zeros",
                "description": "Initial velocity vector q̇(0).",
                "example": [0, 0],
            },
            "t_eval": {
                "type": "float or list[float]",
                "required": False,
                "default": None,
                "description": (
                    "Time(s) at which to evaluate the physical response. "
                    "Results stored as q{i}_at_t{t_val} in final_answer."
                ),
                "example": [0.0, 1.0],
            },
            "compute": {
                "type": "list[str]",
                "required": False,
                "default": "all",
                "description": (
                    "List of computation stages to include. "
                    "Default includes all stages."
                ),
            },
        }
