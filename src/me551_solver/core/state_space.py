"""State-space 변환 솔버 – 1차 상태공간 방정식, 전이행렬, 강제/자유 응답.

Supports two modes:
  1. "state_matrix"  -- direct A, B, x0 input
  2. "second_order"  -- auto-convert Mq̈ + Cq̇ + Kq = Bf(t) to state-space
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _format_matrix(mat: np.ndarray, name: str = "") -> str:
    """Pretty-format a numpy matrix for display."""
    lines = []
    if name:
        lines.append(f"{name} =")
    for row in np.atleast_2d(mat):
        entries = []
        for x in row:
            if np.iscomplex(x) and abs(x.imag) > 1e-12:
                entries.append(f"{x}")
            else:
                entries.append(f"{x.real: .6g}")
        row_str = "  [" + ", ".join(entries) + "]"
        lines.append(row_str)
    return "\n".join(lines)


def _format_vector(vec: np.ndarray, name: str = "") -> str:
    """Pretty-format a numpy vector for display."""
    entries = []
    for x in vec:
        if np.iscomplex(x) and abs(x.imag) > 1e-12:
            entries.append(f"{x}")
        else:
            entries.append(f"{x.real: .6g}")
    vec_str = "[" + ", ".join(entries) + "]"
    if name:
        return f"{name} = {vec_str}"
    return vec_str


def _format_complex(z: complex) -> str:
    """Format a complex number for display."""
    r, i = z.real, z.imag
    if abs(i) < 1e-12:
        return f"{r:.6g}"
    if abs(r) < 1e-12:
        return f"{i:+.6g}j"
    return f"{r:.6g}{i:+.6g}j"


def _format_coeff(val: float) -> str:
    """Format a coefficient for symbolic expressions, using fractions when clean.

    Recognizes patterns like p*sqrt(q)/r for common irrational coefficients.
    """
    frac = _to_fraction(val)
    if frac is not None:
        if frac == "1":
            return ""
        if frac == "-1":
            return "-"
        return frac
    # Try to recognize a*sqrt(b)/c pattern
    sqrt_frac = _to_sqrt_fraction(val)
    if sqrt_frac is not None:
        return sqrt_frac
    return f"{val:.6g}"


def _format_freq(val: float) -> str:
    """Format a frequency value, recognizing sqrt(n) patterns."""
    # Check if val^2 is a clean integer
    val_sq = val * val
    val_sq_round = round(val_sq)
    if abs(val_sq - val_sq_round) < 1e-10 and val_sq_round > 0:
        # Check if val itself is integer
        if abs(val - round(val)) < 1e-10:
            return f"{round(val)}"
        return f"sqrt({val_sq_round})"
    return f"{val:.6g}"


def _to_fraction(val: float, max_denom: int = 1000) -> str | None:
    """Try to represent val as a clean fraction string. Return None if not clean."""
    from fractions import Fraction
    try:
        frac = Fraction(val).limit_denominator(max_denom)
        # Check if the fraction is close enough to the original value
        if abs(float(frac) - val) < 1e-10:
            if frac.denominator == 1:
                return str(frac.numerator)
            return f"({frac.numerator}/{frac.denominator})"
    except (ValueError, OverflowError, ZeroDivisionError):
        pass
    return None


def _to_sqrt_fraction(val: float, max_search: int = 50) -> str | None:
    """Try to express val as (p/q)*sqrt(r) where p,q,r are small integers.

    Returns string like '(2*sqrt(5)/45)' or '-sqrt(3)/2', or None if no match.
    """
    from fractions import Fraction
    sign = 1 if val >= 0 else -1
    absval = abs(val)
    if absval < 1e-15:
        return None

    # Try sqrt(r) for r = 2..max_search
    for r in range(2, max_search + 1):
        sqrt_r = absval * absval / r  # (val/sqrt(r))^2 * r / r ... no
        # val = (p/q) * sqrt(r)  =>  val^2 = (p/q)^2 * r  =>  val^2 / r = (p/q)^2
        ratio_sq = (absval * absval) / r
        ratio = ratio_sq ** 0.5
        frac = Fraction(ratio).limit_denominator(1000)
        if abs(float(frac) - ratio) < 1e-10 and frac.denominator <= 100:
            p = frac.numerator
            q = frac.denominator
            # Verify: (p/q)*sqrt(r) == absval
            check = (p / q) * (r ** 0.5)
            if abs(check - absval) < 1e-10:
                sign_str = "-" if sign < 0 else ""
                if p == 1 and q == 1:
                    return f"{sign_str}sqrt({r})"
                elif q == 1:
                    return f"{sign_str}{p}*sqrt({r})"
                elif p == 1:
                    return f"{sign_str}(sqrt({r})/{q})"
                else:
                    return f"({sign_str}{p}*sqrt({r})/{q})"
    return None


def _format_fraction_or_decimal(val: float) -> str:
    """Format a value as fraction if clean, else decimal."""
    frac = _to_fraction(val)
    if frac is not None:
        return frac
    return f"{val:.6g}"


def _classify_eigenvalue(lam: complex) -> str:
    """Classify a single eigenvalue for stability."""
    r = lam.real
    i = lam.imag
    if abs(r) < 1e-10:
        if abs(i) < 1e-10:
            return "zero (neutral)"
        return "purely imaginary (marginally stable)"
    elif r < 0:
        if abs(i) < 1e-10:
            return "real negative (stable node)"
        return "complex with Re < 0 (stable focus)"
    else:
        if abs(i) < 1e-10:
            return "real positive (unstable node)"
        return "complex with Re > 0 (unstable focus)"


def _classify_system_stability(eigenvalues: np.ndarray) -> dict[str, str]:
    """Classify overall system stability from eigenvalues."""
    real_parts = eigenvalues.real
    max_re = np.max(real_parts)

    if max_re > 1e-10:
        return {
            "stability": "Unstable",
            "classification": "UNSTABLE",
            "details": (
                f"At least one eigenvalue has positive real part "
                f"(max Re(lambda) = {max_re:.6g}). "
                f"Perturbations grow exponentially."
            ),
        }

    # Check for repeated eigenvalues on imaginary axis
    on_axis = eigenvalues[np.abs(real_parts) < 1e-10]
    if len(on_axis) > 0:
        # Round to detect repeated values
        rounded = np.round(on_axis, decimals=8)
        unique, counts = np.unique(rounded, return_counts=True)
        has_repeated = np.any(counts > 1)
        if has_repeated:
            return {
                "stability": "Unstable",
                "classification": "UNSTABLE",
                "details": (
                    "Repeated eigenvalues on the imaginary axis — "
                    "secular (algebraic) growth is possible."
                ),
            }
        return {
            "stability": "Marginally stable",
            "classification": "MARGINALLY STABLE",
            "details": (
                "All eigenvalues have non-positive real part with "
                "purely imaginary eigenvalues being simple. "
                "Bounded oscillations."
            ),
        }

    return {
        "stability": "Asymptotically stable",
        "classification": "STABLE",
        "details": (
            f"All eigenvalues have strictly negative real parts "
            f"(max Re(lambda) = {max_re:.6g}). "
            f"Perturbations decay exponentially to zero."
        ),
    }


def _second_order_to_state_space(
    M: np.ndarray,
    C: np.ndarray,
    K: np.ndarray,
    B_force: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Convert Mq̈ + Cq̇ + Kq = B_force f(t) to state-space form.

    State vector x = [q; q̇].
    Returns (A, B_ss) where:
        A = [[0, I], [-M⁻¹K, -M⁻¹C]]
        B_ss = [[0], [M⁻¹ B_force]]  (or None if B_force is None)
    """
    n = M.shape[0]
    M_inv = np.linalg.inv(M)

    A = np.zeros((2 * n, 2 * n))
    A[:n, n:] = np.eye(n)
    A[n:, :n] = -M_inv @ K
    A[n:, n:] = -M_inv @ C

    B_ss = None
    if B_force is not None:
        B_force = np.atleast_2d(B_force)
        if B_force.shape[0] != n:
            # Try transpose
            if B_force.shape[1] == n:
                B_force = B_force.T
            else:
                raise ValueError(
                    f"B_force shape {B_force.shape} incompatible with {n}-DOF system."
                )
        m_inputs = B_force.shape[1]
        B_ss = np.zeros((2 * n, m_inputs))
        B_ss[n:, :] = M_inv @ B_force

    return A, B_ss


def _compute_residue_matrices(
    eigenvalues: np.ndarray,
    eigvecs: np.ndarray,
) -> list[np.ndarray]:
    """Compute residue matrices R_i = v_i w_i^H where V^{-1} rows are w_i^H.

    For Phi(t) = sum_i R_i exp(lambda_i t).
    """
    try:
        V_inv = np.linalg.inv(eigvecs)
    except np.linalg.LinAlgError:
        return []

    n = len(eigenvalues)
    residues = []
    for i in range(n):
        v_i = eigvecs[:, i].reshape(-1, 1)  # right eigenvector (column)
        w_i_H = V_inv[i, :].reshape(1, -1)  # left eigenvector row (row of V^-1)
        R_i = v_i @ w_i_H
        residues.append(R_i)
    return residues


def _parse_forcing(forcing_desc: Any) -> Callable | None:
    """Parse a forcing function description into a callable f(t) -> array.

    Supported formats:
        - None: no forcing
        - callable: used directly
        - dict with "type": "step", "amplitude": [...], "t_start": 0
        - dict with "type": "harmonic", "amplitude": [...], "omega": w
        - dict with "type": "impulse" (approximated as narrow pulse)
    """
    if forcing_desc is None:
        return None
    if callable(forcing_desc):
        return forcing_desc

    if isinstance(forcing_desc, dict):
        f_type = forcing_desc.get("type", "step")
        amplitude = np.array(forcing_desc.get("amplitude", [1.0]))
        t_start = forcing_desc.get("t_start", 0.0)

        if f_type == "step":
            def f_step(t):
                return amplitude if t >= t_start else np.zeros_like(amplitude)
            return f_step

        elif f_type == "harmonic":
            omega_f = forcing_desc.get("omega", 1.0)
            def f_harmonic(t):
                return amplitude * np.sin(omega_f * t)
            return f_harmonic

        elif f_type == "impulse":
            dt = forcing_desc.get("duration", 0.01)
            def f_impulse(t):
                if t_start <= t <= t_start + dt:
                    return amplitude / dt
                return np.zeros_like(amplitude)
            return f_impulse

        else:
            raise ValueError(f"Unknown forcing type: {f_type}")

    raise ValueError(f"Cannot parse forcing description: {forcing_desc}")


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------


class StateSpaceSolver(BaseSolver):
    """1차 상태공간 방정식의 해석 — 전이행렬, 고유값 분석, 응답 계산.

    Handles both direct state-matrix input and automatic conversion
    from 2nd-order systems (Mq̈ + Cq̇ + Kq = Bf(t)).
    Supports unsymmetric matrices, complex eigenvalues, and forced response.
    """

    def solve(self, params: dict) -> SolverResult:
        """Perform state-space analysis.

        Parameters:
            params: dict — see module docstring for full specification.
                mode: "state_matrix" or "second_order"
        """
        mode = params.get("mode", "state_matrix")
        if mode == "state_matrix":
            return self._solve_state_matrix(params)
        elif mode == "second_order":
            return self._solve_second_order(params)
        else:
            raise ValueError(f"Unknown mode: {mode}")

    # -------------------------------------------------------------------
    # Mode 1: Direct state matrix
    # -------------------------------------------------------------------
    def _solve_state_matrix(self, params: dict) -> SolverResult:
        A_input = params.get("A")
        if A_input is None:
            raise ValueError("'A' (state matrix) is required.")

        A = np.array(A_input, dtype=float)
        n = A.shape[0]

        B_input = params.get("B")
        B = np.array(B_input, dtype=float) if B_input is not None else None

        x0_input = params.get("x0")
        x0 = np.array(x0_input, dtype=float) if x0_input is not None else np.zeros(n)

        t_span = params.get("t_span", [0, 10])
        t_eval_input = params.get("t_eval")
        forcing = params.get("forcing")
        compute = params.get("compute", [
            "eigenvalues", "eigenvectors", "stability",
            "transition_matrix", "response",
        ])

        return self._core_analysis(
            A=A, B=B, x0=x0,
            t_span=t_span, t_eval=t_eval_input,
            forcing=forcing, compute=compute,
            second_order_info=None, params=params,
        )

    # -------------------------------------------------------------------
    # Mode 2: From 2nd-order system
    # -------------------------------------------------------------------
    def _solve_second_order(self, params: dict) -> SolverResult:
        M_input = params.get("M")
        K_input = params.get("K")
        if M_input is None or K_input is None:
            raise ValueError("'M' and 'K' are required for second_order mode.")

        M = np.array(M_input, dtype=float)
        K = np.array(K_input, dtype=float)
        n_dof = M.shape[0]

        C_input = params.get("C")
        C = np.array(C_input, dtype=float) if C_input is not None else np.zeros_like(M)

        B_force_input = params.get("B_force")
        B_force = np.array(B_force_input, dtype=float) if B_force_input is not None else None

        q0_input = params.get("q0")
        q0 = np.array(q0_input, dtype=float) if q0_input is not None else np.zeros(n_dof)

        qdot0_input = params.get("qdot0")
        qdot0 = np.array(qdot0_input, dtype=float) if qdot0_input is not None else np.zeros(n_dof)

        x0 = np.concatenate([q0, qdot0])

        A, B_ss = _second_order_to_state_space(M, C, K, B_force)

        t_span = params.get("t_span", [0, 10])
        t_eval_input = params.get("t_eval")
        forcing = params.get("forcing")
        compute = params.get("compute", [
            "state_matrix", "input_matrix",
            "eigenvalues", "eigenvectors", "stability",
            "transition_matrix", "response",
        ])

        second_order_info = {
            "M": M, "C": C, "K": K,
            "B_force": B_force,
            "n_dof": n_dof,
            "q0": q0, "qdot0": qdot0,
        }

        return self._core_analysis(
            A=A, B=B_ss, x0=x0,
            t_span=t_span, t_eval=t_eval_input,
            forcing=forcing, compute=compute,
            second_order_info=second_order_info, params=params,
        )

    # -------------------------------------------------------------------
    # Core analysis engine
    # -------------------------------------------------------------------
    def _core_analysis(
        self,
        A: np.ndarray,
        B: np.ndarray | None,
        x0: np.ndarray,
        t_span: list,
        t_eval: list | np.ndarray | None,
        forcing: Any,
        compute: list[str],
        second_order_info: dict | None,
        params: dict,
    ) -> SolverResult:
        n = A.shape[0]
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}
        sanity_parts: list[str] = []

        is_second_order = second_order_info is not None

        # ---------------------------------------------------------------
        # 1. Display state matrix A
        # ---------------------------------------------------------------
        if "state_matrix" in compute or True:
            if is_second_order:
                info = second_order_info
                steps.append((
                    "Second-order to state-space conversion",
                    "Given: Mq̈ + Cq̇ + Kq = B_force f(t)\n"
                    "State vector: x = [q; q̇]\n"
                    "State equation: ẋ = Ax + Bu\n"
                    "where A = [[0, I], [-M⁻¹K, -M⁻¹C]],  "
                    "B = [[0], [M⁻¹ B_force]]"
                ))
                steps.append(("Mass matrix M", _format_matrix(info["M"], "M")))
                steps.append(("Damping matrix C", _format_matrix(info["C"], "C")))
                steps.append(("Stiffness matrix K", _format_matrix(info["K"], "K")))

            steps.append((
                "State matrix A",
                _format_matrix(A, "A")
            ))
            final_answer["A"] = A.tolist()

        # ---------------------------------------------------------------
        # 2. Display input matrix B
        # ---------------------------------------------------------------
        if ("input_matrix" in compute or "state_matrix" in compute) and B is not None:
            steps.append((
                "Input matrix B",
                _format_matrix(B, "B")
            ))
            final_answer["B"] = B.tolist()

        # ---------------------------------------------------------------
        # 3. Initial state
        # ---------------------------------------------------------------
        steps.append((
            "Initial state vector x(0)",
            _format_vector(x0, "x(0)")
        ))

        # ---------------------------------------------------------------
        # 4. Eigenvalues of A
        # ---------------------------------------------------------------
        eigenvalues, eigvecs = np.linalg.eig(A)

        # Sort by real part, then imaginary part
        sort_idx = np.lexsort((eigenvalues.imag, eigenvalues.real))
        eigenvalues = eigenvalues[sort_idx]
        eigvecs = eigvecs[:, sort_idx]

        if "eigenvalues" in compute:
            eig_strs = []
            for i, lam in enumerate(eigenvalues):
                eig_type = _classify_eigenvalue(lam)
                eig_strs.append(
                    f"lambda_{i+1} = {_format_complex(lam)}  ({eig_type})"
                )
            steps.append((
                "Eigenvalues of A",
                "\n".join(eig_strs)
            ))
            final_answer["eigenvalues"] = [complex(lam) for lam in eigenvalues]

        # ---------------------------------------------------------------
        # 5. Eigenvectors of A
        # ---------------------------------------------------------------
        if "eigenvectors" in compute:
            evec_strs = []
            for i in range(n):
                evec_strs.append(
                    f"v_{i+1} = " + _format_vector(eigvecs[:, i])
                )
            steps.append((
                "Right eigenvectors of A",
                "\n".join(evec_strs)
            ))
            final_answer["eigenvectors"] = eigvecs.tolist()

        # ---------------------------------------------------------------
        # 6. Stability analysis
        # ---------------------------------------------------------------
        if "stability" in compute:
            stab = _classify_system_stability(eigenvalues)
            steps.append((
                "Stability analysis",
                f"Classification: {stab['stability']}\n{stab['details']}"
            ))
            final_answer["stability"] = stab["stability"]
            final_answer["stability_classification"] = stab["classification"]

        # ---------------------------------------------------------------
        # 7. Transition matrix Phi(t) = e^(At)
        # ---------------------------------------------------------------
        if "transition_matrix" in compute:
            self._compute_transition_matrix(
                A, eigenvalues, eigvecs, n,
                t_eval, steps, final_answer, sanity_parts,
            )

        # ---------------------------------------------------------------
        # 8. Homogeneous / forced response
        # ---------------------------------------------------------------
        if "response" in compute:
            self._compute_response(
                A, B, x0, eigenvalues, eigvecs,
                t_span, t_eval, forcing,
                second_order_info, steps, final_answer, sanity_parts,
            )

        # ---------------------------------------------------------------
        # 9. Physical response reconstruction (for second_order mode)
        # ---------------------------------------------------------------
        if is_second_order and "response" in compute:
            n_dof = second_order_info["n_dof"]
            if "response_x" in final_answer:
                x_resp = np.array(final_answer["response_x"])  # (n_state, n_times)
                q_resp = x_resp[:n_dof, :]
                qdot_resp = x_resp[n_dof:, :]
                final_answer["response_q"] = q_resp.tolist()
                final_answer["response_qdot"] = qdot_resp.tolist()
                steps.append((
                    "Physical response reconstruction",
                    f"Displacement q(t) = x_1..{n_dof}(t)\n"
                    f"Velocity q̇(t) = x_{n_dof+1}..{2*n_dof}(t)"
                ))

        # ---------------------------------------------------------------
        # Sanity checks
        # ---------------------------------------------------------------
        # Phi(0) = I
        Phi_0 = expm(A * 0.0)
        phi0_check = np.allclose(Phi_0, np.eye(n), atol=1e-12)
        sanity_parts.append(
            f"Phi(0) = I: {'PASS' if phi0_check else 'FAIL'}"
        )

        # det(Phi(t)) = e^(tr(A)*t)  — check at t=1
        Phi_1 = expm(A * 1.0)
        det_phi_1 = np.linalg.det(Phi_1)
        expected_det = np.exp(np.trace(A) * 1.0)
        det_check = abs(det_phi_1 - expected_det) < max(1e-6, 1e-6 * abs(expected_det))
        sanity_parts.append(
            f"det(Phi(1)) = {det_phi_1:.6g}, e^(tr(A)*1) = {expected_det:.6g}: "
            f"{'PASS' if det_check else 'FAIL'}"
        )

        # Eigenvalue consistency
        eig_sum = np.sum(eigenvalues)
        trace_A = np.trace(A)
        trace_check = abs(eig_sum - trace_A) < 1e-8
        sanity_parts.append(
            f"Sum of eigenvalues = {_format_complex(eig_sum)}, "
            f"tr(A) = {trace_A:.6g}: {'PASS' if trace_check else 'FAIL'}"
        )

        sanity_check_str = "\n".join(sanity_parts)

        # Build given dict
        given_dict: dict[str, Any] = {"mode": params.get("mode", "state_matrix")}
        if is_second_order:
            given_dict.update({
                "M": second_order_info["M"].tolist(),
                "C": second_order_info["C"].tolist(),
                "K": second_order_info["K"].tolist(),
                "q0": second_order_info["q0"].tolist(),
                "qdot0": second_order_info["qdot0"].tolist(),
                "n_dof": second_order_info["n_dof"],
            })
        else:
            given_dict["A"] = A.tolist()
            given_dict["x0"] = x0.tolist()

        return SolverResult(
            problem_type="State-Space Analysis",
            given=given_dict,
            steps=steps,
            final_answer=final_answer,
            sanity_check=sanity_check_str,
        )

    # -------------------------------------------------------------------
    # Transition matrix computation
    # -------------------------------------------------------------------
    def _compute_transition_matrix(
        self,
        A: np.ndarray,
        eigenvalues: np.ndarray,
        eigvecs: np.ndarray,
        n: int,
        t_eval: list | np.ndarray | None,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
        sanity_parts: list[str],
    ) -> None:
        """Compute and display transition matrix Phi(t) = e^(At)."""

        # --- Eigendecomposition form ---
        is_diagonalizable = True
        try:
            V_inv = np.linalg.inv(eigvecs)
            # Verify: V Lambda V^{-1} ~ A
            Lambda = np.diag(eigenvalues)
            A_reconstructed = eigvecs @ Lambda @ V_inv
            if not np.allclose(A_reconstructed, A, atol=1e-8):
                is_diagonalizable = False
        except np.linalg.LinAlgError:
            is_diagonalizable = False

        if is_diagonalizable:
            steps.append((
                "Transition matrix (eigendecomposition form)",
                "A is diagonalizable: A = V Lambda V⁻¹\n"
                "Phi(t) = V diag(e^(lambda_i * t)) V⁻¹\n"
                "       = sum_i R_i e^(lambda_i * t)\n"
                "where R_i = v_i * w_i^H (residue matrices)"
            ))

            # Compute residue matrices for small systems
            if n <= 6:
                residues = _compute_residue_matrices(eigenvalues, eigvecs)
                if residues:
                    res_strs = []
                    for i, R_i in enumerate(residues):
                        # Show real part if eigenvalue is real
                        if abs(eigenvalues[i].imag) < 1e-12:
                            R_display = R_i.real
                        else:
                            R_display = R_i
                        res_strs.append(
                            f"R_{i+1} (lambda = {_format_complex(eigenvalues[i])}):\n"
                            + _format_matrix(R_display)
                        )
                    steps.append((
                        "Residue matrices R_i",
                        "\n\n".join(res_strs)
                    ))
                    final_answer["residue_matrices"] = [R.tolist() for R in residues]

            # Verify: sum of residues = I
            if n <= 6:
                residues = _compute_residue_matrices(eigenvalues, eigvecs)
                if residues:
                    R_sum = sum(residues)
                    res_sum_check = np.allclose(R_sum, np.eye(n), atol=1e-8)
                    sanity_parts.append(
                        f"Sum of residue matrices = I: "
                        f"{'PASS' if res_sum_check else 'FAIL'}"
                    )
        else:
            steps.append((
                "Transition matrix",
                "A is not diagonalizable (defective matrix).\n"
                "Using scipy.linalg.expm for numeric computation of Phi(t) = e^(At)."
            ))

        # --- Symbolic expression for small systems ---
        if n <= 3 and is_diagonalizable:
            self._symbolic_transition(
                eigenvalues, eigvecs, n, steps, final_answer
            )

        # --- Symbolic Phi(t) via Laplace: L^{-1}{(sI - A)^{-1}} ---
        if n <= 4:
            self._symbolic_transition_via_laplace(
                A, n, steps, final_answer
            )

        # --- Numeric evaluation at specific times ---
        if t_eval is not None:
            t_arr = np.atleast_1d(np.array(t_eval, dtype=float))
            phi_at_t = {}
            for t_val in t_arr:
                Phi_t = expm(A * t_val)
                phi_at_t[f"t={t_val}"] = Phi_t.tolist()
            final_answer["transition_matrix_values"] = phi_at_t
            # Show first and last
            if len(t_arr) <= 5:
                phi_strs = []
                for t_val in t_arr:
                    Phi_t = expm(A * t_val)
                    phi_strs.append(
                        f"Phi({t_val}) =\n" + _format_matrix(Phi_t)
                    )
                steps.append((
                    "Transition matrix at evaluation times",
                    "\n\n".join(phi_strs)
                ))
            else:
                Phi_first = expm(A * t_arr[0])
                Phi_last = expm(A * t_arr[-1])
                steps.append((
                    "Transition matrix at evaluation times",
                    f"Phi({t_arr[0]}) =\n{_format_matrix(Phi_first)}\n\n"
                    f"... ({len(t_arr)} time points) ...\n\n"
                    f"Phi({t_arr[-1]}) =\n{_format_matrix(Phi_last)}"
                ))

    def _symbolic_transition(
        self,
        eigenvalues: np.ndarray,
        eigvecs: np.ndarray,
        n: int,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
    ) -> None:
        """Generate symbolic transition matrix expression for small systems."""
        # Group complex conjugate pairs
        processed = set()
        terms_desc = []

        for i, lam in enumerate(eigenvalues):
            if i in processed:
                continue

            if abs(lam.imag) < 1e-12:
                # Real eigenvalue
                terms_desc.append(
                    f"e^({lam.real:.6g} t) * R_{i+1}"
                )
                processed.add(i)
            else:
                # Complex — find conjugate
                conj_idx = None
                for j in range(i + 1, n):
                    if j not in processed and abs(lam - eigenvalues[j].conj()) < 1e-10:
                        conj_idx = j
                        break

                sigma = lam.real
                omega = abs(lam.imag)
                if conj_idx is not None:
                    terms_desc.append(
                        f"e^({sigma:.6g} t) * "
                        f"[2 Re(R_{i+1}) cos({omega:.6g} t) "
                        f"- 2 Im(R_{i+1}) sin({omega:.6g} t)]"
                    )
                    processed.add(i)
                    processed.add(conj_idx)
                else:
                    terms_desc.append(
                        f"e^({_format_complex(lam)} t) * R_{i+1}"
                    )
                    processed.add(i)

        sym_expr = "Phi(t) = " + "\n       + ".join(terms_desc)
        steps.append((
            "Transition matrix (symbolic form)",
            sym_expr
        ))
        final_answer["transition_matrix_symbolic"] = sym_expr

    # -------------------------------------------------------------------
    # Transition matrix via Laplace transform: Phi(t) = L^{-1}{(sI-A)^{-1}}
    # -------------------------------------------------------------------
    def _symbolic_transition_via_laplace(
        self,
        A: np.ndarray,
        n: int,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
    ) -> None:
        """Compute Phi(t) = L^{-1}{(sI - A)^{-1}} symbolically via sympy."""
        try:
            import sympy as sp
        except ImportError:
            return

        s = sp.Symbol("s")
        t_sym = sp.Symbol("t", positive=True)

        # Build symbolic A
        A_sym = sp.Matrix([[sp.Rational(x).limit_denominator(10000)
                            for x in row] for row in A.tolist()])

        sI = s * sp.eye(n)
        sI_minus_A = sI - A_sym

        # Step 1: Show (sI - A)
        si_strs = []
        for i in range(n):
            row_entries = [str(sI_minus_A[i, j]) for j in range(n)]
            si_strs.append("  [" + ",  ".join(row_entries) + "]")
        step1_text = "(sI - A) =\n" + "\n".join(si_strs)

        # Step 2: Characteristic polynomial = det(sI - A)
        char_poly = sp.expand(sI_minus_A.det())
        step2_text = f"det(sI - A) = {char_poly}"

        # Step 3: Adjugate matrix
        adj = sI_minus_A.adjugate()
        adj_strs = []
        for i in range(n):
            row_entries = [str(sp.expand(adj[i, j])) for j in range(n)]
            adj_strs.append("  [" + ",  ".join(row_entries) + "]")
        step3_text = "adj(sI - A) =\n" + "\n".join(adj_strs)

        # Step 4: (sI - A)^{-1} entries as rational functions
        resolvent_strs = []
        for i in range(n):
            for j in range(n):
                entry = sp.simplify(adj[i, j] / char_poly)
                resolvent_strs.append(f"  [{i+1},{j+1}]: {entry}")
        step4_text = "(sI - A)⁻¹ = adj(sI - A) / det(sI - A):\n" + "\n".join(resolvent_strs)

        steps.append((
            "Transition matrix via Laplace: Φ(t) = L⁻¹{(sI - A)⁻¹}",
            step1_text + "\n\n" + step2_text + "\n\n" + step3_text + "\n\n" + step4_text
        ))

        # Step 5: Inverse Laplace of each entry
        Phi_sym = sp.zeros(n, n)
        ilt_strs = []
        for i in range(n):
            for j in range(n):
                F_s = adj[i, j] / char_poly
                f_t = sp.inverse_laplace_transform(F_s, s, t_sym)
                # Clean up Heaviside
                f_t = f_t.rewrite(sp.exp)
                f_t = f_t.replace(sp.Heaviside(t_sym), sp.Integer(1))
                f_t = f_t.replace(sp.Heaviside, lambda *args: sp.Integer(1))
                f_t = sp.simplify(sp.trigsimp(f_t.rewrite(sp.cos)))
                Phi_sym[i, j] = f_t
                ilt_strs.append(f"  Φ_{i+1}{j+1}(t) = {f_t}")

        steps.append((
            "Inverse Laplace of each entry → Φ(t)",
            "\n".join(ilt_strs)
        ))

        # Store symbolic Phi(t) as formatted matrix
        phi_matrix_strs = []
        for i in range(n):
            row = [str(Phi_sym[i, j]) for j in range(n)]
            phi_matrix_strs.append("  [" + ",  ".join(row) + "]")
        phi_matrix_text = "Φ(t) =\n" + "\n".join(phi_matrix_strs)

        steps.append((
            "Complete transition matrix Φ(t) (via Laplace)",
            phi_matrix_text
        ))
        final_answer["transition_matrix_laplace"] = phi_matrix_text

    # -------------------------------------------------------------------
    # Response computation
    # -------------------------------------------------------------------
    def _compute_response(
        self,
        A: np.ndarray,
        B: np.ndarray | None,
        x0: np.ndarray,
        eigenvalues: np.ndarray,
        eigvecs: np.ndarray,
        t_span: list,
        t_eval: list | np.ndarray | None,
        forcing: Any,
        second_order_info: dict | None,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
        sanity_parts: list[str],
    ) -> None:
        """Compute free and/or forced response."""
        n = A.shape[0]

        # Build time array
        if t_eval is not None:
            t_arr = np.atleast_1d(np.array(t_eval, dtype=float))
        else:
            t_start, t_end = float(t_span[0]), float(t_span[1])
            n_pts = max(200, int(50 * (t_end - t_start)))
            n_pts = min(n_pts, 2000)
            t_arr = np.linspace(t_start, t_end, n_pts)

        f_func = _parse_forcing(forcing)

        if f_func is None:
            # --- Homogeneous response: x(t) = Phi(t) x(0) ---
            steps.append((
                "Homogeneous response",
                "x(t) = Phi(t) x(0) = e^(At) x(0)\n"
                "No forcing function provided."
            ))

            # Compute response at all time points
            x_response = np.zeros((n, len(t_arr)))
            for k, t_val in enumerate(t_arr):
                x_response[:, k] = expm(A * t_val) @ x0

            # Symbolic expression via eigendecomposition (for small systems)
            if n <= 4:
                self._symbolic_free_response(
                    eigenvalues, eigvecs, x0, n, steps, final_answer
                )

        else:
            # --- Forced response via solve_ivp ---
            if B is None:
                raise ValueError(
                    "Input matrix B is required for forced response computation."
                )

            steps.append((
                "Forced response",
                "x(t) = Phi(t) x(0) + integral_0^t Phi(t-tau) B f(tau) d_tau\n"
                "Computed numerically via scipy.integrate.solve_ivp (RK45)."
            ))

            # --- Analytical closed-form for step / impulse inputs ---
            if n <= 6 and forcing is not None and isinstance(forcing, dict):
                self._symbolic_forced_response(
                    A, B, x0, eigenvalues, eigvecs, n, forcing,
                    steps, final_answer,
                )

            def ode_rhs(t, x):
                f_val = np.atleast_1d(f_func(t)).flatten()
                return A @ x + B @ f_val

            sol = solve_ivp(
                ode_rhs,
                [t_arr[0], t_arr[-1]],
                x0,
                t_eval=t_arr,
                method="RK45",
                rtol=1e-10,
                atol=1e-12,
                max_step=0.01 * (t_arr[-1] - t_arr[0]),
            )

            if not sol.success:
                steps.append((
                    "Integration warning",
                    f"solve_ivp message: {sol.message}"
                ))

            x_response = sol.y  # shape (n, len(t_arr))
            t_arr = sol.t

        # Store results
        final_answer["response_t"] = t_arr.tolist()
        final_answer["response_x"] = x_response.tolist()

        # Summary of response at key time points
        resp_summary_strs = []
        sample_indices = [0]
        if len(t_arr) > 1:
            sample_indices.append(len(t_arr) // 4)
            sample_indices.append(len(t_arr) // 2)
            sample_indices.append(3 * len(t_arr) // 4)
            sample_indices.append(len(t_arr) - 1)
        sample_indices = sorted(set(sample_indices))

        for idx in sample_indices:
            t_val = t_arr[idx]
            x_val = x_response[:, idx]
            resp_summary_strs.append(
                f"x(t={t_val:.4g}) = {_format_vector(x_val)}"
            )

        steps.append((
            "Response summary at selected times",
            "\n".join(resp_summary_strs)
        ))

        # Verify initial condition
        x_at_0 = x_response[:, 0]
        ic_check = np.allclose(x_at_0, x0, atol=1e-6)
        sanity_parts.append(
            f"Response x(0) matches initial condition: "
            f"{'PASS' if ic_check else 'FAIL'}"
        )

    def _symbolic_free_response(
        self,
        eigenvalues: np.ndarray,
        eigvecs: np.ndarray,
        x0: np.ndarray,
        n: int,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
    ) -> None:
        """Generate symbolic free response expression for small systems."""
        try:
            V_inv = np.linalg.inv(eigvecs)
        except np.linalg.LinAlgError:
            return

        # Modal initial conditions: alpha = V^{-1} x0
        alpha = V_inv @ x0
        term_strs = []

        processed = set()
        for i in range(n):
            if i in processed:
                continue
            lam = eigenvalues[i]
            a = alpha[i]

            if abs(lam.imag) < 1e-12:
                # Real eigenvalue
                lam_r = lam.real
                a_r = a.real
                if abs(a_r) > 1e-14:
                    term_strs.append(
                        f"{a_r:.6g} * v_{i+1} * e^({lam_r:.6g} t)"
                    )
                processed.add(i)
            else:
                # Complex pair
                conj_idx = None
                for j in range(i + 1, n):
                    if j not in processed and abs(lam - eigenvalues[j].conj()) < 1e-10:
                        conj_idx = j
                        break

                sigma = lam.real
                omega = abs(lam.imag)

                if conj_idx is not None:
                    # Use eigenvalue with positive imaginary part
                    # 2 Re(alpha_+ v_+ e^(lambda_+ t))
                    # = 2 e^(sigma t) [Re(alpha_+ v_+) cos(omega t) - Im(alpha_+ v_+) sin(omega t)]
                    if lam.imag > 0:
                        av = a * eigvecs[:, i]
                    else:
                        av = alpha[conj_idx] * eigvecs[:, conj_idx]
                    re_av = av.real
                    im_av = av.imag
                    term_strs.append(
                        f"2 e^({sigma:.6g} t) * "
                        f"[{_format_vector(re_av)} cos({omega:.6g} t) "
                        f"- {_format_vector(im_av)} sin({omega:.6g} t)]"
                    )
                    processed.add(i)
                    processed.add(conj_idx)
                else:
                    if abs(a) > 1e-14:
                        term_strs.append(
                            f"{_format_complex(a)} * v_{i+1} * "
                            f"e^({_format_complex(lam)} t)"
                        )
                    processed.add(i)

        if term_strs:
            sym_resp = "x(t) = " + "\n     + ".join(term_strs)
        else:
            sym_resp = "x(t) = 0 (zero initial conditions)"

        steps.append((
            "Free response (modal expansion)",
            sym_resp
        ))
        final_answer["free_response_symbolic"] = sym_resp

    # -------------------------------------------------------------------
    # Analytical forced response (closed-form via transition matrix)
    # -------------------------------------------------------------------
    def _symbolic_forced_response(
        self,
        A: np.ndarray,
        B: np.ndarray,
        x0: np.ndarray,
        eigenvalues: np.ndarray,
        eigvecs: np.ndarray,
        n: int,
        forcing: dict,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
    ) -> None:
        """Derive closed-form forced response using the transition matrix.

        For unit step u(t):
          x(t) = Phi(t) x(0) + A^{-1}(Phi(t) - I) B
               = Phi(t) x(0) + sum_i (R_i B / lambda_i)(e^{lambda_i t} - 1)
          where steady-state = -A^{-1} B  (for stable systems)

        For impulse delta(t):
          x(t) = Phi(t) x(0) + Phi(t) B   (for t > 0)
               = Phi(t) (x(0) + B)
        """
        f_type = forcing.get("type", "")
        F0 = forcing.get("F0", 1.0)

        if f_type not in ("step", "impulse"):
            return

        try:
            V_inv = np.linalg.inv(eigvecs)
        except np.linalg.LinAlgError:
            return

        residues = _compute_residue_matrices(eigenvalues, eigvecs)
        if not residues:
            return

        # Effective forcing amplitude vector through B
        # B is (n, m), for scalar input m=1
        B_col = B[:, 0] * F0 if B.shape[1] >= 1 else B.flatten() * F0

        if f_type == "impulse":
            # x(t) = Phi(t) (x(0) + B*F0) for t > 0
            x0_eff = x0 + B_col
            steps.append((
                "Analytical impulse response",
                f"For f(t) = {F0}*delta(t):\n"
                f"  x(t) = Phi(t) * [x(0) + B*F0]   (t > 0)\n"
                f"  Effective IC: x(0+) = {_format_vector(x0_eff)}\n"
                "  (Reduces to free response with modified initial condition)"
            ))
            # Delegate to symbolic free response with modified IC
            self._symbolic_free_response(
                eigenvalues, eigvecs, x0_eff, n, steps, final_answer
            )
            return

        # --- Step response: f(t) = F0 * u(t) ---
        # Convolution: integral_0^t Phi(t-tau) B F0 dtau
        #   = sum_i R_i B F0 * (1/lambda_i)(e^{lambda_i t} - 1)
        #
        # Total: x(t) = Phi(t) x(0) + sum_i (R_i B F0 / lambda_i)(e^{lambda_i t} - 1)

        # Steady-state value (if A is invertible and system is stable)
        try:
            A_inv = np.linalg.inv(A)
            x_ss = -A_inv @ B_col
        except np.linalg.LinAlgError:
            x_ss = None

        # Build symbolic expression for each state variable
        # Coefficient of constant term: sum_i (-R_i B / lambda_i) = -A^{-1} B = x_ss
        # Coefficient of e^{lambda_i t}: (R_i B / lambda_i) + contribution from Phi(t) x(0)
        #   From free: alpha_i * v_i  where alpha = V^{-1} x(0)
        #   From forced: R_i B / lambda_i

        alpha = V_inv @ x0  # modal ICs for free response

        # For each eigenvalue, compute the total coefficient of e^{lambda_i t}
        # c_i = alpha_i * v_i + (R_i @ B_col) / lambda_i
        coeffs = []
        for i in range(n):
            lam = eigenvalues[i]
            free_part = alpha[i] * eigvecs[:, i]  # from Phi(t) x(0)
            if abs(lam) > 1e-14:
                forced_part = (residues[i] @ B_col) / lam
            else:
                forced_part = np.zeros(n, dtype=complex)
            coeffs.append(free_part + forced_part)

        # Build human-readable expression
        derivation_lines = []
        derivation_lines.append(
            f"f(t) = {F0} * u(t)  (unit step scaled by {F0})"
        )
        derivation_lines.append("")
        derivation_lines.append(
            "Convolution integral:")
        derivation_lines.append(
            "  integral_0^t Phi(t-tau) B f(tau) dtau"
        )
        derivation_lines.append(
            "  = sum_i (R_i B / lambda_i) (e^(lambda_i t) - 1)"
        )
        derivation_lines.append("")

        if x_ss is not None:
            derivation_lines.append(
                f"Steady-state: x_ss = -A^(-1) B * {F0} = {_format_vector(x_ss)}"
            )
            derivation_lines.append("")

        # Build per-component expressions
        derivation_lines.append("Combining free + forced response for each state:")
        derivation_lines.append("")

        # Variable names for display
        var_names = [f"x_{k+1}(t)" for k in range(n)]

        processed = set()
        for k in range(n):
            terms = []
            # Constant term (steady-state contribution)
            if x_ss is not None:
                const_val = x_ss[k].real
                if abs(const_val) > 1e-14:
                    terms.append(_format_fraction_or_decimal(const_val))

            # Exponential terms grouped by conjugate pairs
            processed_inner = set()
            for i in range(n):
                if i in processed_inner:
                    continue
                lam = eigenvalues[i]
                c_i = coeffs[i][k]  # coefficient for state k

                if abs(lam.imag) < 1e-12:
                    # Real eigenvalue
                    c_r = c_i.real
                    lam_r = lam.real
                    if abs(c_r) > 1e-14:
                        terms.append(
                            f"{_format_coeff(c_r)} e^({lam_r:.6g} t)"
                        )
                    processed_inner.add(i)
                else:
                    # Complex conjugate pair
                    conj_idx = None
                    for j in range(i + 1, n):
                        if j not in processed_inner and abs(lam - eigenvalues[j].conj()) < 1e-10:
                            conj_idx = j
                            break

                    sigma = lam.real
                    omega = abs(lam.imag)

                    if conj_idx is not None:
                        # Use the eigenvalue with POSITIVE imaginary part
                        # so that: 2 Re(c_+ e^{lambda_+ t})
                        # = 2 e^{sigma t} [Re(c_+) cos(omega t) - Im(c_+) sin(omega t)]
                        if lam.imag > 0:
                            c_pos = c_i
                        else:
                            c_pos = coeffs[conj_idx][k]
                        cos_coeff = 2.0 * c_pos.real
                        sin_coeff = -2.0 * c_pos.imag

                        if abs(cos_coeff) > 1e-14:
                            terms.append(
                                f"{_format_coeff(cos_coeff)} e^({sigma:.6g} t) cos({_format_freq(omega)} t)"
                            )
                        if abs(sin_coeff) > 1e-14:
                            terms.append(
                                f"{_format_coeff(sin_coeff)} e^({sigma:.6g} t) sin({_format_freq(omega)} t)"
                            )
                        processed_inner.add(i)
                        processed_inner.add(conj_idx)
                    else:
                        if abs(c_i) > 1e-14:
                            terms.append(
                                f"{_format_complex(c_i)} e^({_format_complex(lam)} t)"
                            )
                        processed_inner.add(i)

            if terms:
                expr = " + ".join(terms).replace("+ -", "- ")
                derivation_lines.append(f"  {var_names[k]} = {expr}")
            else:
                derivation_lines.append(f"  {var_names[k]} = 0")

        steps.append((
            "Analytical step response (closed-form)",
            "\n".join(derivation_lines)
        ))
        final_answer["forced_response_symbolic"] = "\n".join(derivation_lines)

    # -------------------------------------------------------------------
    # Input template
    # -------------------------------------------------------------------
    def get_input_template(self) -> dict:
        """필요한 입력 파라미터 설명을 반환한다."""
        return {
            "mode": {
                "type": "str",
                "required": False,
                "default": "state_matrix",
                "description": (
                    "Analysis mode: 'state_matrix' (direct A, B, x0) "
                    "or 'second_order' (auto-convert from M, C, K)."
                ),
            },
            "A": {
                "type": "list[list[float]]",
                "required": "for state_matrix mode",
                "description": "State matrix A (n x n).",
                "example": [[0, 1], [-4, -2]],
            },
            "B": {
                "type": "list[list[float]]",
                "required": False,
                "description": "Input matrix B (n x m). Required if forcing is provided.",
                "example": [[0], [1]],
            },
            "x0": {
                "type": "list[float]",
                "required": False,
                "default": "zeros",
                "description": "Initial state vector.",
                "example": [1, 0],
            },
            "M": {
                "type": "list[list[float]]",
                "required": "for second_order mode",
                "description": "Mass matrix (n x n).",
                "example": [[1, 0], [0, 1]],
            },
            "C": {
                "type": "list[list[float]]",
                "required": False,
                "default": "zeros",
                "description": "Damping matrix (n x n).",
                "example": [[2, -1], [-1, 2]],
            },
            "K": {
                "type": "list[list[float]]",
                "required": "for second_order mode",
                "description": "Stiffness matrix (n x n).",
                "example": [[4, -2], [-2, 4]],
            },
            "B_force": {
                "type": "list[list[float]]",
                "required": False,
                "description": "Force input matrix for second_order mode.",
                "example": [[1], [0]],
            },
            "q0": {
                "type": "list[float]",
                "required": False,
                "default": "zeros",
                "description": "Initial displacement (second_order mode).",
                "example": [1, 0],
            },
            "qdot0": {
                "type": "list[float]",
                "required": False,
                "default": "zeros",
                "description": "Initial velocity (second_order mode).",
                "example": [0, 0],
            },
            "t_span": {
                "type": "[float, float]",
                "required": False,
                "default": [0, 10],
                "description": "Time range for response computation.",
                "example": [0, 10],
            },
            "t_eval": {
                "type": "list[float]",
                "required": False,
                "default": None,
                "description": (
                    "Specific times at which to evaluate the transition matrix "
                    "and response."
                ),
                "example": [0, 0.1, 0.2, 0.5, 1.0],
            },
            "forcing": {
                "type": "dict or callable",
                "required": False,
                "default": None,
                "description": (
                    "Forcing function description. Supported types:\n"
                    "  - None: free response\n"
                    "  - callable f(t) -> array of size m\n"
                    "  - dict with 'type': 'step'/'harmonic'/'impulse' and 'amplitude'"
                ),
                "example": {"type": "step", "amplitude": [1.0], "t_start": 0},
            },
            "compute": {
                "type": "list[str]",
                "required": False,
                "default": "all",
                "description": (
                    "List of computations: 'state_matrix', 'input_matrix', "
                    "'eigenvalues', 'eigenvectors', 'stability', "
                    "'transition_matrix', 'response'."
                ),
            },
        }
