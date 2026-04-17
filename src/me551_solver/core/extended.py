"""Extended solver for gyroscopic and nonconservative systems.

Handles the general equation of motion:
    M q̈ + (C + G) q̇ + (K + H) q = f(t)

where:
    M = mass matrix (symmetric positive definite)
    K = stiffness matrix (symmetric)
    C = damping matrix (symmetric, optional)
    G = gyroscopic matrix (skew-symmetric, G^T = -G, optional)
    H = nonconservative/circulatory matrix (skew-symmetric, H^T = -H, optional)

Provides:
    - Left and right eigenvector computation
    - Bi-orthogonality verification
    - Complex mode analysis
    - Gyroscopic stabilization detection
    - Flutter / divergence instability classification
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.linalg import inv, eig, norm

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _format_matrix(mat: np.ndarray, name: str = "") -> str:
    """Pretty-format a numpy matrix (real or complex) for display."""
    lines = []
    if name:
        lines.append(f"{name} =")
    for row in mat:
        entries = []
        for x in row:
            if np.iscomplex(x) and abs(x.imag) > 1e-12:
                entries.append(f"{x.real: .6g}{x.imag:+.6g}j")
            else:
                entries.append(f"{np.real(x): .6g}")
        row_str = "  [" + ", ".join(entries) + "]"
        lines.append(row_str)
    return "\n".join(lines)


def _format_vector(vec: np.ndarray, name: str = "") -> str:
    """Pretty-format a numpy vector (real or complex) for display."""
    entries = []
    for x in vec:
        if np.iscomplex(x) and abs(x.imag) > 1e-12:
            entries.append(f"{x.real: .6g}{x.imag:+.6g}j")
        else:
            entries.append(f"{np.real(x): .6g}")
    vec_str = "[" + ", ".join(entries) + "]"
    if name:
        return f"{name} = {vec_str}"
    return vec_str


def _format_complex(z: complex) -> str:
    """Format a single complex number."""
    if abs(z.imag) < 1e-12:
        return f"{z.real:.6g}"
    if abs(z.real) < 1e-12:
        return f"{z.imag:+.6g}j"
    return f"{z.real:.6g}{z.imag:+.6g}j"


def _is_skew_symmetric(A: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if A is skew-symmetric: A^T = -A."""
    return np.allclose(A, -A.T, atol=tol)


def _is_symmetric(A: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if A is symmetric: A^T = A."""
    return np.allclose(A, A.T, atol=tol)


def _sort_eigenvalues_by_frequency(eigenvalues: np.ndarray,
                                    right_vecs: np.ndarray,
                                    left_vecs: np.ndarray):
    """Sort eigenvalues by imaginary part (ascending frequency).

    Returns sorted eigenvalues, right eigenvectors, left eigenvectors.
    """
    idx = np.argsort(np.abs(eigenvalues.imag))
    return eigenvalues[idx], right_vecs[:, idx], left_vecs[:, idx]


def _pair_conjugate_eigenvalues(eigenvalues: np.ndarray):
    """Group eigenvalues into conjugate pairs.

    Returns list of tuples:
        For conjugate pair: (index_pos_imag, index_neg_imag)
        For real eigenvalue: (index,)
    """
    n = len(eigenvalues)
    used = set()
    pairs = []

    for i in range(n):
        if i in used:
            continue
        ev = eigenvalues[i]
        if abs(ev.imag) < 1e-10:
            # Real eigenvalue
            pairs.append((i,))
            used.add(i)
        else:
            # Find conjugate partner
            best_j = -1
            best_err = np.inf
            for j in range(n):
                if j in used or j == i:
                    continue
                err = abs(ev - eigenvalues[j].conjugate())
                if err < best_err:
                    best_err = err
                    best_j = j
            if best_j >= 0 and best_err < 1e-8:
                # Ensure the one with positive imaginary part comes first
                if ev.imag > 0:
                    pairs.append((i, best_j))
                else:
                    pairs.append((best_j, i))
                used.add(i)
                used.add(best_j)
            else:
                # Unpaired complex eigenvalue (should not happen for real A)
                pairs.append((i,))
                used.add(i)
    return pairs


def _classify_system(G: np.ndarray | None, H: np.ndarray | None,
                     C: np.ndarray | None) -> str:
    """Classify the system type based on which matrices are present and nonzero."""
    has_G = G is not None and norm(G) > 1e-14
    has_H = H is not None and norm(H) > 1e-14
    has_C = C is not None and norm(C) > 1e-14

    if has_G and has_H:
        return "mixed_gyroscopic_nonconservative"
    elif has_G and not has_H:
        if has_C:
            return "damped_gyroscopic"
        return "pure_gyroscopic"
    elif has_H and not has_G:
        if has_C:
            return "damped_nonconservative"
        return "pure_nonconservative"
    elif has_C:
        return "damped_only"
    else:
        return "conservative"


def _classify_stability_from_eigenvalues(eigenvalues: np.ndarray) -> dict[str, Any]:
    """Classify stability from complex eigenvalues of the state matrix.

    Returns dict with:
        stability: str -- overall stability verdict
        stable_modes: int -- count of stable modes
        unstable_modes: int -- count of unstable modes
        flutter_modes: list[int] -- indices of flutter-unstable modes
        divergence_modes: list[int] -- indices of divergence-unstable modes
        gyroscopic_stable: bool -- True if all eigenvalues are purely imaginary
        details: str
    """
    tol = 1e-10
    real_parts = eigenvalues.real
    imag_parts = eigenvalues.imag
    max_re = np.max(real_parts)

    flutter_modes = []
    divergence_modes = []
    stable_count = 0

    for i, ev in enumerate(eigenvalues):
        if ev.real > tol:
            if abs(ev.imag) > tol:
                flutter_modes.append(i)
            else:
                divergence_modes.append(i)
        else:
            stable_count += 1

    # Check if all eigenvalues are purely imaginary (gyroscopic stable)
    all_imaginary = np.all(np.abs(real_parts) < tol)

    if max_re > tol:
        if flutter_modes and not divergence_modes:
            stability = "Flutter unstable"
            details = (
                f"Flutter instability detected: {len(flutter_modes)} eigenvalue(s) "
                f"with positive real part and nonzero imaginary part.\n"
                f"Max Re(lambda) = {max_re:.6g}.\n"
                f"The system oscillates with exponentially growing amplitude."
            )
        elif divergence_modes and not flutter_modes:
            stability = "Divergence unstable"
            details = (
                f"Divergence instability detected: {len(divergence_modes)} real "
                f"positive eigenvalue(s).\n"
                f"Max Re(lambda) = {max_re:.6g}.\n"
                f"The system diverges monotonically from equilibrium."
            )
        else:
            stability = "Unstable (flutter + divergence)"
            details = (
                f"Both flutter ({len(flutter_modes)} modes) and divergence "
                f"({len(divergence_modes)} modes) instabilities present.\n"
                f"Max Re(lambda) = {max_re:.6g}."
            )
    elif all_imaginary:
        stability = "Gyroscopically stable"
        details = (
            "All eigenvalues are purely imaginary (on the imaginary axis).\n"
            "The system exhibits bounded oscillations without energy dissipation.\n"
            "This is characteristic of a gyroscopic conservative system."
        )
    elif max_re < -tol:
        stability = "Asymptotically stable"
        details = (
            f"All eigenvalues have strictly negative real parts.\n"
            f"Max Re(lambda) = {max_re:.6g}.\n"
            f"Perturbations decay exponentially to zero."
        )
    else:
        # Eigenvalues on imaginary axis but not all purely imaginary
        stability = "Marginally stable"
        details = (
            "All eigenvalues have non-positive real parts with "
            "purely imaginary eigenvalues being simple.\n"
            "Bounded oscillations."
        )

    return {
        "stability": stability,
        "stable_modes": stable_count,
        "unstable_modes": len(flutter_modes) + len(divergence_modes),
        "flutter_modes": flutter_modes,
        "divergence_modes": divergence_modes,
        "gyroscopic_stable": all_imaginary,
        "max_real_part": float(max_re),
        "details": details,
    }


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

class ExtendedSolver(BaseSolver):
    """Solver for gyroscopic and nonconservative vibration systems.

    Handles the general equation of motion:
        M q'' + (C + G) q' + (K + H) q = 0

    where G is skew-symmetric (gyroscopic) and H is skew-symmetric
    (nonconservative / circulatory).

    Computes state-space formulation, left/right eigenvectors,
    bi-orthogonality, complex modes, and stability classification.
    """

    def solve(self, params: dict) -> SolverResult:
        """Perform extended eigenanalysis for gyroscopic/nonconservative systems.

        Parameters:
            params: dict with keys:
                mode: "gyroscopic", "nonconservative", or "general_unsymmetric"
                M: mass matrix (symmetric positive definite)
                K: stiffness matrix (symmetric)
                G: gyroscopic matrix (skew-symmetric, optional)
                H: nonconservative/circulatory matrix (skew-symmetric, optional)
                C: damping matrix (symmetric, optional)
                initial_q: initial displacement (optional)
                initial_qdot: initial velocity (optional)
                compute: list of stages to include (optional, default all)
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

        G_input = params.get("G")
        H_input = params.get("H")
        C_input = params.get("C")

        G = np.array(G_input, dtype=float) if G_input is not None else np.zeros((n_dof, n_dof))
        H = np.array(H_input, dtype=float) if H_input is not None else np.zeros((n_dof, n_dof))
        C = np.array(C_input, dtype=float) if C_input is not None else np.zeros((n_dof, n_dof))

        initial_q = params.get("initial_q")
        initial_qdot = params.get("initial_qdot")
        q0 = np.array(initial_q, dtype=float) if initial_q is not None else np.zeros(n_dof)
        qdot0 = np.array(initial_qdot, dtype=float) if initial_qdot is not None else np.zeros(n_dof)

        compute = params.get("compute", [
            "state_matrix", "eigenvalues", "right_eigenvectors",
            "left_eigenvectors", "bi_orthogonality", "complex_modes",
            "stability", "response",
        ])

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # ---------------------------------------------------------------
        # 1. System identification -- validate matrix properties
        # ---------------------------------------------------------------
        checks = []
        m_sym = _is_symmetric(M)
        k_sym = _is_symmetric(K)
        g_skew = _is_skew_symmetric(G)
        h_skew = _is_skew_symmetric(H)
        c_sym = _is_symmetric(C)

        checks.append(f"M symmetric: {'PASS' if m_sym else 'FAIL'}")
        checks.append(f"K symmetric: {'PASS' if k_sym else 'FAIL'}")
        if G_input is not None:
            checks.append(f"G skew-symmetric (G^T = -G): {'PASS' if g_skew else 'FAIL'}")
        if H_input is not None:
            checks.append(f"H skew-symmetric (H^T = -H): {'PASS' if h_skew else 'FAIL'}")
        if C_input is not None:
            checks.append(f"C symmetric: {'PASS' if c_sym else 'FAIL'}")

        # Check M positive definite
        m_eigvals = np.linalg.eigvalsh(M)
        m_pd = np.all(m_eigvals > 0)
        checks.append(f"M positive definite: {'PASS' if m_pd else 'FAIL'} "
                       f"(min eigenvalue = {np.min(m_eigvals):.6g})")

        sys_type = _classify_system(G_input and G, H_input and H, C_input and C)
        # Reclassify using actual norm checks
        sys_type = _classify_system(G, H, C)

        sys_type_display = {
            "pure_gyroscopic": "Pure Gyroscopic (M q'' + G q' + K q = 0)",
            "damped_gyroscopic": "Damped Gyroscopic (M q'' + (C+G) q' + K q = 0)",
            "pure_nonconservative": "Nonconservative (M q'' + (K+H) q = 0)",
            "damped_nonconservative": "Damped Nonconservative (M q'' + C q' + (K+H) q = 0)",
            "mixed_gyroscopic_nonconservative": "Mixed (M q'' + (C+G) q' + (K+H) q = 0)",
            "damped_only": "Damped (M q'' + C q' + K q = 0)",
            "conservative": "Conservative (M q'' + K q = 0)",
        }.get(sys_type, sys_type)

        checks.append(f"\nSystem classification: {sys_type_display}")
        final_answer["system_type"] = sys_type

        steps.append(("System identification and matrix property checks",
                       "\n".join(checks)))

        # Display all matrices
        mat_strs = [_format_matrix(M, "M"), _format_matrix(K, "K")]
        if norm(G) > 1e-14:
            mat_strs.append(_format_matrix(G, "G"))
        if norm(H) > 1e-14:
            mat_strs.append(_format_matrix(H, "H"))
        if norm(C) > 1e-14:
            mat_strs.append(_format_matrix(C, "C"))
        steps.append(("System matrices", "\n\n".join(mat_strs)))

        # ---------------------------------------------------------------
        # 2. State-space formulation
        # ---------------------------------------------------------------
        if "state_matrix" in compute:
            M_inv = inv(M)

            # Effective velocity matrix: C + G
            V_eff = C + G
            # Effective stiffness matrix: K + H
            K_eff = K + H

            # State matrix A: x = [q, q_dot]^T, x_dot = A x
            # A = [[    0,        I     ],
            #      [-M^{-1}(K+H), -M^{-1}(C+G)]]
            A = np.zeros((2 * n_dof, 2 * n_dof))
            A[:n_dof, n_dof:] = np.eye(n_dof)
            A[n_dof:, :n_dof] = -M_inv @ K_eff
            A[n_dof:, n_dof:] = -M_inv @ V_eff

            steps.append((
                "State-space formulation: x = [q, q_dot]^T, x_dot = A x",
                f"A = [[0, I], [-M^(-1)(K+H), -M^(-1)(C+G)]]\n\n"
                + _format_matrix(A, "A")
            ))
            final_answer["state_matrix"] = A.tolist()
        else:
            # Still need A for further computations
            M_inv = inv(M)
            V_eff = C + G
            K_eff = K + H
            A = np.zeros((2 * n_dof, 2 * n_dof))
            A[:n_dof, n_dof:] = np.eye(n_dof)
            A[n_dof:, :n_dof] = -M_inv @ K_eff
            A[n_dof:, n_dof:] = -M_inv @ V_eff

        # ---------------------------------------------------------------
        # 3. Right eigenvectors: A phi_r = lambda_r phi_r
        # ---------------------------------------------------------------
        if "eigenvalues" in compute or "right_eigenvectors" in compute:
            eigvals_right, phi_right = eig(A)

            # Sort by imaginary part (ascending frequency)
            sort_idx = np.argsort(np.abs(eigvals_right.imag))
            eigvals_right = eigvals_right[sort_idx]
            phi_right = phi_right[:, sort_idx]
        else:
            eigvals_right, phi_right = eig(A)
            sort_idx = np.argsort(np.abs(eigvals_right.imag))
            eigvals_right = eigvals_right[sort_idx]
            phi_right = phi_right[:, sort_idx]

        if "eigenvalues" in compute:
            eig_strs = []
            for i, ev in enumerate(eigvals_right):
                eig_strs.append(f"lambda_{i+1} = {_format_complex(ev)}")
            steps.append(("Eigenvalues of state matrix A", "\n".join(eig_strs)))
            final_answer["eigenvalues"] = [complex(ev) for ev in eigvals_right]

        if "right_eigenvectors" in compute:
            rev_strs = []
            for i in range(2 * n_dof):
                rev_strs.append(
                    f"phi_R,{i+1} = {_format_vector(phi_right[:, i])}"
                )
            steps.append((
                "Right eigenvectors: A phi_R = lambda phi_R",
                "\n".join(rev_strs)
            ))
            final_answer["right_eigenvectors"] = phi_right.tolist()

        # ---------------------------------------------------------------
        # 4. Left eigenvectors: psi^T A = lambda psi^T  =>  A^T psi = lambda psi
        # ---------------------------------------------------------------
        if "left_eigenvectors" in compute:
            eigvals_left, psi_left_raw = eig(A.T)

            # Match left eigenvectors to right eigenvectors by eigenvalue
            # The eigenvalues of A^T are the same as A, but possibly in
            # different order. We reorder left eigenvectors to match the
            # right eigenvalue ordering.
            psi_left = np.zeros_like(psi_left_raw)
            used_left = set()
            for i, ev_r in enumerate(eigvals_right):
                best_j = -1
                best_err = np.inf
                for j in range(2 * n_dof):
                    if j in used_left:
                        continue
                    err = abs(eigvals_left[j] - ev_r)
                    if err < best_err:
                        best_err = err
                        best_j = j
                if best_j >= 0:
                    psi_left[:, i] = psi_left_raw[:, best_j]
                    used_left.add(best_j)

            # Normalize left eigenvectors so that psi_i^T phi_i = 1
            for i in range(2 * n_dof):
                dot_val = psi_left[:, i] @ phi_right[:, i]
                if abs(dot_val) > 1e-14:
                    psi_left[:, i] = psi_left[:, i] / dot_val

            lev_strs = []
            for i in range(2 * n_dof):
                lev_strs.append(
                    f"psi_L,{i+1} = {_format_vector(psi_left[:, i])}"
                )
            steps.append((
                "Left eigenvectors: psi^T A = lambda psi^T  (eigenvectors of A^T)",
                "Normalized so that psi_i^T phi_i = 1\n\n" + "\n".join(lev_strs)
            ))
            final_answer["left_eigenvectors"] = psi_left.tolist()

        # ---------------------------------------------------------------
        # 5. Bi-orthogonality verification: psi_i^T phi_j = delta_ij
        # ---------------------------------------------------------------
        if "bi_orthogonality" in compute and "left_eigenvectors" in compute:
            bi_orth_matrix = psi_left.T @ phi_right

            # Check: should be close to identity
            identity_check = np.allclose(bi_orth_matrix, np.eye(2 * n_dof), atol=1e-6)

            bi_orth_detail = _format_matrix(bi_orth_matrix, "psi^T * phi")
            bi_orth_detail += (
                f"\n\nBi-orthogonality check (should be identity): "
                f"{'PASS' if identity_check else 'FAIL'}"
            )

            # Also check psi_i^T A phi_j = lambda_j delta_ij
            psiT_A_phi = psi_left.T @ A @ phi_right
            spectral_diag = np.diag(eigvals_right)
            spectral_check = np.allclose(psiT_A_phi, spectral_diag, atol=1e-6)
            bi_orth_detail += (
                f"\n\npsi^T A phi = diag(lambda) check: "
                f"{'PASS' if spectral_check else 'FAIL'}"
            )

            steps.append(("Bi-orthogonality verification", bi_orth_detail))
            final_answer["bi_orthogonality_matrix"] = bi_orth_matrix.tolist()
            final_answer["bi_orthogonality_pass"] = bool(identity_check)

        # ---------------------------------------------------------------
        # 6. Complex mode analysis
        # ---------------------------------------------------------------
        if "complex_modes" in compute:
            pairs = _pair_conjugate_eigenvalues(eigvals_right)

            mode_strs = []
            mode_data = []
            mode_num = 0

            for pair in pairs:
                if len(pair) == 2:
                    # Conjugate pair -- one complex mode
                    i_pos, i_neg = pair
                    ev = eigvals_right[i_pos]
                    mode_num += 1

                    sigma = -ev.real       # decay rate
                    omega_d = abs(ev.imag)  # damped natural frequency
                    omega_n = abs(ev)       # undamped natural frequency magnitude
                    if omega_n > 1e-14:
                        zeta = sigma / omega_n
                    else:
                        zeta = 0.0

                    # Mode shape: upper half of right eigenvector (displacement part)
                    mode_shape = phi_right[:n_dof, i_pos]

                    mode_info = {
                        "mode_number": mode_num,
                        "eigenvalue": complex(ev),
                        "omega_n": float(omega_n),
                        "omega_d": float(omega_d),
                        "zeta": float(zeta),
                        "sigma": float(sigma),
                        "mode_shape_real": mode_shape.real.tolist(),
                        "mode_shape_imag": mode_shape.imag.tolist(),
                    }
                    mode_data.append(mode_info)

                    mode_str = (
                        f"--- Complex Mode {mode_num} ---\n"
                        f"  Eigenvalue: {_format_complex(ev)} (conjugate: {_format_complex(ev.conjugate())})\n"
                        f"  omega_n = |lambda| = {omega_n:.6g} rad/s\n"
                        f"  omega_d = |Im(lambda)| = {omega_d:.6g} rad/s\n"
                        f"  zeta = -Re(lambda)/|lambda| = {zeta:.6g}\n"
                        f"  sigma (decay rate) = {sigma:.6g} 1/s\n"
                        f"  Mode shape (displacement):\n"
                        f"    Real part: {_format_vector(mode_shape.real)}\n"
                        f"    Imag part: {_format_vector(mode_shape.imag)}"
                    )
                    mode_strs.append(mode_str)
                else:
                    # Real eigenvalue
                    i_real = pair[0]
                    ev = eigvals_right[i_real]
                    if abs(ev.imag) > 1e-10:
                        # Unpaired complex -- treat as a mode
                        mode_num += 1
                        omega_n = abs(ev)
                        mode_shape = phi_right[:n_dof, i_real]
                        mode_info = {
                            "mode_number": mode_num,
                            "eigenvalue": complex(ev),
                            "omega_n": float(omega_n),
                            "omega_d": float(abs(ev.imag)),
                            "zeta": float(-ev.real / omega_n) if omega_n > 1e-14 else 0.0,
                            "sigma": float(-ev.real),
                            "mode_shape_real": mode_shape.real.tolist(),
                            "mode_shape_imag": mode_shape.imag.tolist(),
                        }
                        mode_data.append(mode_info)
                    else:
                        # Pure real eigenvalue
                        mode_num += 1
                        mode_shape = phi_right[:n_dof, i_real]
                        mode_info = {
                            "mode_number": mode_num,
                            "eigenvalue": complex(ev),
                            "omega_n": float(abs(ev.real)),
                            "omega_d": 0.0,
                            "zeta": 1.0 if ev.real < 0 else -1.0,
                            "sigma": float(-ev.real),
                            "mode_shape_real": mode_shape.real.tolist(),
                            "mode_shape_imag": mode_shape.imag.tolist(),
                        }
                        mode_data.append(mode_info)
                        status = "decaying" if ev.real < 0 else "divergent"
                        mode_str = (
                            f"--- Real Mode {mode_num} ---\n"
                            f"  Eigenvalue: {_format_complex(ev)} ({status})\n"
                            f"  Mode shape (displacement):\n"
                            f"    Real part: {_format_vector(mode_shape.real)}"
                        )
                        mode_strs.append(mode_str)

            steps.append((
                "Complex mode analysis",
                "\n\n".join(mode_strs) if mode_strs else "No modes extracted."
            ))
            final_answer["complex_modes"] = mode_data

            # Collect omega_n and zeta for each mode
            for md in mode_data:
                mn = md["mode_number"]
                final_answer[f"omega_n_{mn}"] = md["omega_n"]
                final_answer[f"zeta_{mn}"] = md["zeta"]

        # ---------------------------------------------------------------
        # 7. Gyroscopic-specific analysis
        # ---------------------------------------------------------------
        if "gyroscopic" in sys_type or sys_type == "conservative":
            gyro_strs = []

            # Energy: G does no work (q_dot^T G q_dot = 0)
            gyro_strs.append(
                "Gyroscopic energy property:\n"
                "  P_gyro = q_dot^T G q_dot = 0  (since G is skew-symmetric)\n"
                "  Gyroscopic forces redistribute energy among modes but do not add or remove energy."
            )

            # Check if eigenvalues are purely imaginary (gyroscopic stability)
            all_imag = np.all(np.abs(eigvals_right.real) < 1e-10)
            if all_imag and norm(G) > 1e-14:
                gyro_strs.append(
                    "\nGyroscopic stability:\n"
                    "  All eigenvalues are purely imaginary.\n"
                    "  The system is gyroscopically stable with bounded oscillations."
                )

            # Gyroscopic stabilization check: test stability without G
            if norm(G) > 1e-14:
                A_noG = np.zeros((2 * n_dof, 2 * n_dof))
                A_noG[:n_dof, n_dof:] = np.eye(n_dof)
                A_noG[n_dof:, :n_dof] = -M_inv @ (K + H)
                A_noG[n_dof:, n_dof:] = -M_inv @ C
                eigvals_noG = np.linalg.eigvals(A_noG)
                max_re_noG = np.max(eigvals_noG.real)

                max_re_withG = np.max(eigvals_right.real)

                if max_re_noG > 1e-10 and max_re_withG < 1e-10:
                    gyro_strs.append(
                        "\nGyroscopic stabilization detected:\n"
                        f"  Without G: max Re(lambda) = {max_re_noG:.6g} (UNSTABLE)\n"
                        f"  With G:    max Re(lambda) = {max_re_withG:.6g} (STABLE)\n"
                        "  The gyroscopic matrix G stabilizes an otherwise unstable system.\n"
                        "  Note: Adding damping can destroy gyroscopic stabilization\n"
                        "  (Kelvin-Tait-Chetaev theorem)."
                    )
                    final_answer["gyroscopic_stabilization"] = True
                elif max_re_noG > 1e-10 and max_re_withG > 1e-10:
                    gyro_strs.append(
                        "\nGyroscopic stabilization insufficient:\n"
                        f"  Without G: max Re(lambda) = {max_re_noG:.6g} (UNSTABLE)\n"
                        f"  With G:    max Re(lambda) = {max_re_withG:.6g} (STILL UNSTABLE)\n"
                        "  Gyroscopic parameter may be below the critical threshold."
                    )
                    final_answer["gyroscopic_stabilization"] = False
                else:
                    final_answer["gyroscopic_stabilization"] = False

            if gyro_strs:
                steps.append(("Gyroscopic-specific analysis", "\n".join(gyro_strs)))

        # ---------------------------------------------------------------
        # 8. Nonconservative-specific analysis
        # ---------------------------------------------------------------
        if "nonconservative" in sys_type and norm(H) > 1e-14:
            nc_strs = []

            # Check if K alone is positive definite
            k_eigvals = np.linalg.eigvalsh(K)
            k_pd = np.all(k_eigvals > 0)

            nc_strs.append(
                f"Stiffness matrix K analysis:\n"
                f"  Eigenvalues of K: {np.array2string(k_eigvals, precision=6)}\n"
                f"  K positive definite: {'YES' if k_pd else 'NO'}"
            )

            # Effective stiffness K_eff = K + H
            keff_eigvals = np.linalg.eigvals(K_eff)
            nc_strs.append(
                f"\nEffective stiffness K+H analysis:\n"
                f"  K+H is generally non-symmetric (due to H^T = -H)\n"
                f"  Eigenvalues of K+H: "
                + ", ".join(_format_complex(ev) for ev in keff_eigvals)
            )

            # Flutter instability mechanism
            has_flutter = any(
                eigvals_right[i].real > 1e-10 and abs(eigvals_right[i].imag) > 1e-10
                for i in range(len(eigvals_right))
            )
            if has_flutter and k_pd:
                nc_strs.append(
                    "\nFlutter instability from nonconservative forces:\n"
                    "  Despite K being positive definite, the circulatory matrix H\n"
                    "  introduces non-symmetric stiffness that causes flutter.\n"
                    "  This is a classic nonconservative instability mechanism."
                )
            elif has_flutter:
                nc_strs.append(
                    "\nFlutter instability detected.\n"
                    "  Nonconservative forces (H) contribute to the instability."
                )

            steps.append(("Nonconservative-specific analysis", "\n".join(nc_strs)))

        # ---------------------------------------------------------------
        # 9. Stability classification
        # ---------------------------------------------------------------
        if "stability" in compute:
            stab = _classify_stability_from_eigenvalues(eigvals_right)
            final_answer["stability"] = stab["stability"]
            final_answer["stability_details"] = stab["details"]
            final_answer["max_real_part"] = stab["max_real_part"]
            final_answer["flutter_mode_count"] = len(stab["flutter_modes"])
            final_answer["divergence_mode_count"] = len(stab["divergence_modes"])
            final_answer["gyroscopic_stable"] = stab["gyroscopic_stable"]

            stab_str = (
                f"Stability verdict: {stab['stability']}\n\n"
                f"{stab['details']}\n\n"
                f"Max Re(lambda) = {stab['max_real_part']:.6g}\n"
                f"Flutter modes: {len(stab['flutter_modes'])}\n"
                f"Divergence modes: {len(stab['divergence_modes'])}"
            )

            # Physical interpretation
            if stab["gyroscopic_stable"]:
                stab_str += (
                    "\n\nPhysical interpretation:\n"
                    "  Purely imaginary eigenvalues indicate constant-amplitude oscillations.\n"
                    "  Energy is conserved and redistributed among modes by gyroscopic coupling."
                )
            elif stab["stability"] == "Asymptotically stable":
                stab_str += (
                    "\n\nPhysical interpretation:\n"
                    "  All perturbations decay exponentially.\n"
                    "  Damping dissipates energy faster than any destabilizing mechanism."
                )
            elif "Flutter" in stab["stability"]:
                stab_str += (
                    "\n\nPhysical interpretation:\n"
                    "  Flutter: two modes coalesce in frequency and exchange stability.\n"
                    "  One mode gains energy from the nonconservative field while oscillating."
                )
            elif "Divergence" in stab["stability"]:
                stab_str += (
                    "\n\nPhysical interpretation:\n"
                    "  Divergence: the system monotonically moves away from equilibrium.\n"
                    "  The effective stiffness is negative in at least one direction."
                )

            steps.append(("Stability classification", stab_str))

        # ---------------------------------------------------------------
        # 10. Free response (if initial conditions given)
        # ---------------------------------------------------------------
        if "response" in compute and (norm(q0) > 1e-14 or norm(qdot0) > 1e-14):
            x0 = np.concatenate([q0, qdot0])

            # Modal decomposition via bi-orthogonality: eta_i(0) = psi_i^T x0
            if "left_eigenvectors" in compute:
                eta0 = psi_left.T @ x0

                resp_strs = []
                resp_strs.append(
                    "Free response via modal expansion:\n"
                    "  x(t) = sum_i eta_i(0) * exp(lambda_i * t) * phi_R,i\n"
                    "  where eta_i(0) = psi_L,i^T * x0"
                )

                ic_strs = []
                for i in range(2 * n_dof):
                    ic_strs.append(
                        f"  eta_{i+1}(0) = {_format_complex(eta0[i])}"
                    )
                resp_strs.append("\nModal initial conditions:\n" + "\n".join(ic_strs))

                # Physical response expression (symbolic form)
                resp_strs.append(
                    "\nPhysical displacement response:\n"
                    "  q_i(t) = sum_r eta_r(0) * exp(lambda_r * t) * phi_R,i,r\n"
                    "  (sum over all 2N state-space modes)"
                )

                # Evaluate at t=0 for verification
                q_check = np.zeros(n_dof, dtype=complex)
                for i in range(2 * n_dof):
                    q_check += eta0[i] * phi_right[:n_dof, i]
                q_check_real = q_check.real

                ic_check = np.allclose(q_check_real, q0, atol=1e-8)
                resp_strs.append(
                    f"\nIC reconstruction check: q(0) from modal expansion = "
                    f"{_format_vector(q_check_real)}\n"
                    f"Given q(0) = {_format_vector(q0)}\n"
                    f"Match: {'PASS' if ic_check else 'FAIL'}"
                )

                steps.append(("Free response analysis", "\n".join(resp_strs)))
                final_answer["modal_ic"] = [complex(e) for e in eta0]

                # Evaluate at specific times if requested
                t_eval = params.get("t_eval")
                if t_eval is not None:
                    t_arr = np.atleast_1d(np.array(t_eval, dtype=float))
                    for t_val in t_arr:
                        x_t = np.zeros(2 * n_dof, dtype=complex)
                        for i in range(2 * n_dof):
                            x_t += eta0[i] * np.exp(eigvals_right[i] * t_val) * phi_right[:, i]
                        q_t = x_t[:n_dof].real
                        for i in range(n_dof):
                            final_answer[f"q{i+1}_at_t{t_val}"] = float(q_t[i])
            else:
                steps.append((
                    "Free response",
                    "Left eigenvectors not computed. Include 'left_eigenvectors' "
                    "in compute list for modal response via bi-orthogonality."
                ))

        # ---------------------------------------------------------------
        # Sanity check
        # ---------------------------------------------------------------
        sanity_parts: list[str] = []

        if m_sym and m_pd:
            sanity_parts.append("M is symmetric positive definite: OK")
        else:
            sanity_parts.append("WARNING: M may not be symmetric positive definite")

        if G_input is not None and not g_skew:
            sanity_parts.append("WARNING: G is NOT skew-symmetric -- check input")
        if H_input is not None and not h_skew:
            sanity_parts.append("WARNING: H is NOT skew-symmetric -- check input")

        sanity_parts.append(f"System type: {sys_type_display}")

        # Eigenvalue count check
        n_eigs = len(eigvals_right)
        sanity_parts.append(
            f"State-space dimension: 2N = {2 * n_dof}, eigenvalue count: {n_eigs}"
        )

        # Stability
        if "stability" in final_answer:
            sanity_parts.append(f"Stability: {final_answer['stability']}")

        # Bi-orthogonality
        if "bi_orthogonality_pass" in final_answer:
            sanity_parts.append(
                f"Bi-orthogonality: {'PASS' if final_answer['bi_orthogonality_pass'] else 'FAIL'}"
            )

        sanity_check_str = "\n".join(sanity_parts)

        return SolverResult(
            problem_type="Extended Analysis (Gyroscopic/Nonconservative)",
            given={
                "M": M.tolist(),
                "K": K.tolist(),
                "G": G.tolist() if G_input is not None else None,
                "H": H.tolist() if H_input is not None else None,
                "C": C.tolist() if C_input is not None else None,
                "n_dof": n_dof,
                "mode": params.get("mode", "general_unsymmetric"),
                "initial_q": q0.tolist(),
                "initial_qdot": qdot0.tolist(),
            },
            steps=steps,
            final_answer=final_answer,
            sanity_check=sanity_check_str,
        )

    def get_input_template(self) -> dict:
        """Return description of required input parameters."""
        return {
            "mode": {
                "type": "str",
                "required": False,
                "default": "general_unsymmetric",
                "description": (
                    "System mode: 'gyroscopic', 'nonconservative', or 'general_unsymmetric'."
                ),
            },
            "M": {
                "type": "list[list[float]] or np.ndarray",
                "required": True,
                "description": "Mass matrix (n x n, symmetric positive definite).",
                "example": [[1, 0], [0, 1]],
            },
            "K": {
                "type": "list[list[float]] or np.ndarray",
                "required": True,
                "description": "Stiffness matrix (n x n, symmetric).",
                "example": [[3, -1], [-1, 3]],
            },
            "G": {
                "type": "list[list[float]] or np.ndarray",
                "required": False,
                "default": "zeros",
                "description": "Gyroscopic matrix (n x n, skew-symmetric: G^T = -G).",
                "example": [[0, 2], [-2, 0]],
            },
            "H": {
                "type": "list[list[float]] or np.ndarray",
                "required": False,
                "default": "zeros",
                "description": (
                    "Nonconservative/circulatory matrix (n x n, skew-symmetric: H^T = -H)."
                ),
                "example": [[0, 1], [-1, 0]],
            },
            "C": {
                "type": "list[list[float]] or np.ndarray",
                "required": False,
                "default": "zeros",
                "description": "Damping matrix (n x n, symmetric).",
                "example": [[0.1, 0], [0, 0.1]],
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
                "description": "Initial velocity vector q_dot(0).",
            },
            "compute": {
                "type": "list[str]",
                "required": False,
                "default": "all",
                "description": (
                    "List of computation stages: 'state_matrix', 'eigenvalues', "
                    "'right_eigenvectors', 'left_eigenvectors', 'bi_orthogonality', "
                    "'complex_modes', 'stability', 'response'."
                ),
            },
            "t_eval": {
                "type": "float or list[float]",
                "required": False,
                "default": None,
                "description": "Time(s) at which to evaluate the response.",
            },
        }
