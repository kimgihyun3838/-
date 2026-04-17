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

        elif f_type == "exponential":
            F0_val = forcing_desc.get("F0", 1.0)
            a_val = forcing_desc.get("a", 1.0)
            def f_exponential(t):
                return amplitude * F0_val * np.exp(-a_val * t)
            return f_exponential

        elif f_type == "exp_truncated":
            F0_val = forcing_desc.get("F0", 1.0)
            a_val = forcing_desc.get("a", 1.0)
            t1_val = forcing_desc.get("t1", 1.0)
            def f_exp_trunc(t):
                if 0 <= t <= t1_val:
                    return amplitude * F0_val * np.exp(-a_val * t)
                return np.zeros_like(amplitude)
            return f_exp_trunc

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
        elif mode == "symbolic":
            return self._solve_symbolic(params)
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
    # Mode 3: Fully symbolic (sympy matrices, no numeric conversion)
    # -------------------------------------------------------------------
    def _solve_symbolic(self, params: dict) -> SolverResult:
        """Fully symbolic state-space analysis — ωn, ζ, ωd 등 문자 그대로 출력.

        params:
            A_sym: sympy.Matrix  (state matrix with symbols)
            B_sym: sympy.Matrix  (input matrix, optional)
            x0_sym: sympy.Matrix (initial state, optional)
            forcing: dict        (same format as numeric mode)
            use_omega_d: bool    (True → introduce ωd = ωn√(1-ζ²) for clean output)
        """
        import sympy as sp
        from sympy import (
            Symbol, Matrix, eye, exp, cos, sin, sqrt, Rational,
            simplify, trigsimp, expand, expand_trig,
            inverse_laplace_transform, Heaviside, nsimplify,
        )

        A_sym: sp.Matrix = params["A_sym"]
        B_sym: sp.Matrix | None = params.get("B_sym")
        x0_sym: sp.Matrix | None = params.get("x0_sym")
        forcing = params.get("forcing")
        use_omega_d = params.get("use_omega_d", True)

        n = A_sym.shape[0]
        if x0_sym is None:
            x0_sym = sp.zeros(n, 1)

        t_s = Symbol("t", positive=True)
        tau_s = Symbol("tau", positive=True)
        s = Symbol("s")

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {"mode": "symbolic"}
        sanity_parts: list[str] = []

        # ── optional ωd substitution for cleaner output ──
        omega_n = Symbol("omega_n", positive=True)
        zeta = Symbol("zeta", positive=True)
        omega_d = Symbol("omega_d", positive=True)
        sigma_sym = Symbol("sigma", positive=True)

        A_work = A_sym
        wd_subs_info = ""
        if use_omega_d:
            # Replace omega_n^2 → sigma^2 + omega_d^2, 2*zeta*omega_n → 2*sigma
            # where sigma = zeta*omega_n, omega_d = omega_n*sqrt(1-zeta^2)
            if A_sym.has(omega_n) and A_sym.has(zeta):
                A_work = A_sym.subs([
                    (zeta * omega_n, sigma_sym),
                    (omega_n**2, sigma_sym**2 + omega_d**2),
                ])
                # Double check: also try to catch 2*omega_n*zeta patterns
                A_work = A_work.subs(omega_n * zeta, sigma_sym)
                wd_subs_info = (
                    "치환: σ = ζωₙ,  ωd = ωₙ√(1-ζ²),  ωₙ² = σ² + ωd²\n"
                    "  (underdamped 가정: ζ < 1)\n"
                )

        # ── 1. Display matrices ──
        steps.append((
            "State equation",
            f"ẋ = Ax + Bu\n\n"
            f"A = {A_sym}\n"
            + (f"B = {B_sym}\n" if B_sym is not None else "")
            + f"x(0) = {x0_sym.T}\n"
            + (f"\n{wd_subs_info}" if wd_subs_info else "")
            + (f"Working A (치환 후) = {A_work}" if A_work != A_sym else ""),
        ))

        # ── 2. Characteristic polynomial ──
        sI_A = s * eye(n) - A_work
        char_poly = expand(sI_A.det())
        steps.append((
            "Characteristic polynomial",
            f"det(sI - A) = {char_poly}",
        ))
        final_answer["char_poly"] = str(char_poly)

        # ── 3. Eigenvalues ──
        eigenvals = A_work.eigenvals()
        eig_strs = []
        for ev, mult in eigenvals.items():
            ev_simplified = simplify(ev)
            label = f"  λ = {ev_simplified}"
            if mult > 1:
                label += f"  (multiplicity {mult})"
            eig_strs.append(label)
        steps.append(("Eigenvalues", "\n".join(eig_strs)))
        final_answer["eigenvalues"] = [str(ev) for ev in eigenvals.keys()]

        # ── 4. Transition matrix Φ(t) via Laplace ──
        adj = sI_A.adjugate()
        Phi = sp.zeros(n, n)
        phi_strs = []
        for i in range(n):
            for j in range(n):
                F_s = adj[i, j] / char_poly
                f_t = inverse_laplace_transform(F_s, s, t_s)
                f_t = f_t.rewrite(sp.exp)
                f_t = f_t.replace(Heaviside(t_s), sp.Integer(1))
                f_t = f_t.replace(Heaviside, lambda *a: sp.Integer(1))
                f_t = simplify(trigsimp(f_t.rewrite(cos)))
                Phi[i, j] = f_t
                phi_strs.append(f"  Φ_{i + 1}{j + 1}(t) = {f_t}")

        steps.append(("Transition matrix Φ(t) = e^(At)", "\n".join(phi_strs)))
        final_answer["Phi"] = str(Phi)

        # ── 5. Free response (if non-zero IC) ──
        if x0_sym != sp.zeros(n, 1):
            x_free = simplify(Phi * x0_sym)
            free_strs = [
                f"  x_{i + 1}(t) = {x_free[i, 0]}" for i in range(n)
            ]
            steps.append((
                "Free response: x(t) = Φ(t)·x(0)",
                "\n".join(free_strs),
            ))

        # ── 6. Forced response ──
        if forcing and B_sym is not None:
            # Adapt B to same substitution
            B_work = B_sym
            if use_omega_d and B_sym.has(omega_n) and A_sym.has(zeta):
                B_work = B_sym.subs([
                    (zeta * omega_n, sigma_sym),
                    (omega_n**2, sigma_sym**2 + omega_d**2),
                ])

            self._symbolic_forced_full(
                Phi, A_work, B_work, x0_sym, n, t_s, tau_s,
                forcing, steps, final_answer,
            )

        # ── Sanity ──
        # Phi(0) = I check (symbolic)
        Phi_0 = Phi.subs(t_s, 0)
        phi0_ok = simplify(Phi_0 - eye(n)) == sp.zeros(n, n)
        sanity_parts.append(f"Φ(0) = I: {'PASS' if phi0_ok else 'CHECK — ' + str(Phi_0)}")

        return SolverResult(
            problem_type="State-Space Analysis (Symbolic)",
            given={"A": str(A_sym),
                   "B": str(B_sym) if B_sym else "None",
                   "x0": str(x0_sym.T)},
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    def _symbolic_forced_full(
        self,
        Phi: "sp.Matrix",
        A_work: "sp.Matrix",
        B_work: "sp.Matrix",
        x0_sym: "sp.Matrix",
        n: int,
        t_s: "sp.Symbol",
        tau_s: "sp.Symbol",
        forcing: dict,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
    ) -> None:
        """Symbolic forced response via transition matrix convolution.

        Supports: step, exponential, exp_truncated, piecewise (N regions).
        For piecewise, forcing["pieces"] is a list of dicts:
            [{"expr": "F0*exp(-a*tau)", "start": 0, "end": "t1"},
             {"expr": "F1*sin(omega*tau)", "start": "t1", "end": "t2"},
             {"expr": "0", "start": "t2", "end": "inf"}]
        """
        import sympy as sp
        from sympy import simplify, nsimplify, oo

        f_type = forcing.get("type", "")

        def _pp(val, name, **kw):
            if isinstance(val, str):
                return sp.Symbol(name, **kw)
            return nsimplify(val)

        Phi_shifted = Phi.subs(t_s, t_s - tau_s)
        PhiB = Phi_shifted * B_work

        # ── Piecewise forcing (general N-region) ──
        if f_type == "piecewise":
            self._symbolic_piecewise_convolution(
                Phi, PhiB, B_work, x0_sym, n, t_s, tau_s,
                forcing, steps, final_answer,
            )
            return

        # ── Single-expression forcing types ──
        F0 = _pp(forcing.get("F0", 1), "F_0", positive=True)
        a_p = _pp(forcing.get("a", 1), "a", positive=True)

        if f_type == "exponential":
            f_tau = F0 * sp.exp(-a_p * tau_s)
            integrand = PhiB * f_tau

            steps.append((
                "Convolution integral (exponential forcing)",
                f"f(t) = {F0}·e^(-{a_p}·t),  t ≥ 0\n\n"
                f"x(t) = Φ(t)·x₀ + ∫₀ᵗ Φ(t-τ)·B·f(τ) dτ",
            ))

            x_forced = sp.zeros(n, 1)
            r_lines = []
            for i in range(n):
                result = self._sym_integrate_conv(
                    integrand[i, 0], tau_s, sp.S.Zero, t_s,
                )
                x_forced[i, 0] = result

            x_free = Phi * x0_sym
            for i in range(n):
                total = simplify(x_free[i, 0] + x_forced[i, 0])
                r_lines.append(f"  x_{i + 1}(t) = {total}")

            steps.append(("Closed-form response", "\n".join(r_lines)))
            final_answer["forced_response_symbolic"] = "\n".join(r_lines)

        elif f_type == "exp_truncated":
            # Convert to piecewise internally
            t1 = _pp(forcing.get("t1", 1), "t_1", positive=True)
            pw_forcing = {
                "type": "piecewise",
                "pieces": [
                    {"expr": F0 * sp.exp(-a_p * tau_s), "start": sp.S.Zero, "end": t1},
                    {"expr": sp.S.Zero, "start": t1, "end": oo},
                ],
            }
            self._symbolic_piecewise_convolution(
                Phi, PhiB, B_work, x0_sym, n, t_s, tau_s,
                pw_forcing, steps, final_answer,
            )

        elif f_type == "step":
            f_tau = F0
            integrand = PhiB * f_tau
            steps.append((
                "Step response (symbolic)",
                f"f(t) = {F0}·u(t)\n\n"
                f"x(t) = Φ(t)·x₀ + ∫₀ᵗ Φ(t-τ)·B·{F0} dτ",
            ))
            x_forced = sp.zeros(n, 1)
            r_lines = []
            for i in range(n):
                result = self._sym_integrate_conv(
                    integrand[i, 0], tau_s, sp.S.Zero, t_s,
                )
                x_forced[i, 0] = result
            x_free = Phi * x0_sym
            for i in range(n):
                total = simplify(x_free[i, 0] + x_forced[i, 0])
                r_lines.append(f"  x_{i + 1}(t) = {total}")
            steps.append(("Closed-form step response", "\n".join(r_lines)))
            final_answer["forced_response_symbolic"] = "\n".join(r_lines)

    # -------------------------------------------------------------------
    # General piecewise convolution via transition matrix
    # -------------------------------------------------------------------
    def _symbolic_piecewise_convolution(
        self,
        Phi: "sp.Matrix",
        PhiB: "sp.Matrix",
        B_work: "sp.Matrix",
        x0_sym: "sp.Matrix",
        n: int,
        t_s: "sp.Symbol",
        tau_s: "sp.Symbol",
        forcing: dict,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
    ) -> None:
        """N-region piecewise forcing via transition matrix.

        Strategy:
          For each time region k (t_{k-1} < t < t_k):
            x(t) = Φ(t - t_{k-1}) · x(t_{k-1})           (propagate prior state)
                 + ∫_{t_{k-1}}^{t} Φ(t-τ) · B · f_k(τ) dτ  (current forcing)

          At each boundary t_k:
            x(t_k) = evaluated from above

          For regions where f(t) = 0:
            x(t) = Φ(t - t_{k-1}) · x(t_{k-1})           (pure free vibration)
        """
        import sympy as sp
        from sympy import simplify, oo

        pieces = forcing["pieces"]
        # Each piece: {"expr": sympy_expr_in_tau, "start": sym, "end": sym}
        # If expr/start/end are strings, parse them
        _sym_ns = {
            "tau": tau_s, "t": tau_s,
            "F0": sp.Symbol("F_0", positive=True),
            "F_0": sp.Symbol("F_0", positive=True),
            "F1": sp.Symbol("F_1", positive=True),
            "F_1": sp.Symbol("F_1", positive=True),
            "a": sp.Symbol("a", positive=True),
            "b": sp.Symbol("b", positive=True),
            "omega": sp.Symbol("omega", positive=True),
            "w": sp.Symbol("omega", positive=True),
            "t1": sp.Symbol("t_1", positive=True),
            "t_1": sp.Symbol("t_1", positive=True),
            "t2": sp.Symbol("t_2", positive=True),
            "t_2": sp.Symbol("t_2", positive=True),
            "t3": sp.Symbol("t_3", positive=True),
            "t_3": sp.Symbol("t_3", positive=True),
            "m": sp.Symbol("m", positive=True),
            "k": sp.Symbol("k", positive=True),
            "exp": sp.exp, "sin": sp.sin, "cos": sp.cos,
            "sqrt": sp.sqrt, "pi": sp.pi,
        }

        def _to_sym(val):
            if isinstance(val, sp.Basic):
                return val
            s = str(val).strip()
            if s.lower() in ("inf", "oo", "infinity"):
                return oo
            return sp.sympify(s, locals=_sym_ns)

        parsed_pieces: list[tuple[sp.Expr, sp.Expr, sp.Expr]] = []
        desc_lines = []
        for p in pieces:
            expr = _to_sym(p["expr"])
            lo = _to_sym(p["start"])
            hi = _to_sym(p["end"])
            parsed_pieces.append((expr, lo, hi))
            if hi == oo:
                desc_lines.append(f"  f(t) = {expr},  t > {lo}")
            else:
                desc_lines.append(f"  f(t) = {expr},  {lo} < t < {hi}")

        # Collect boundary points (exclude 0 and ∞)
        boundaries: list[sp.Expr] = []
        for _, lo, hi in parsed_pieces:
            if lo != 0 and lo != sp.S.Zero and lo not in boundaries:
                boundaries.append(lo)
            if hi != oo and not hi.is_infinite and hi not in boundaries:
                boundaries.append(hi)
        # Sort by symbol name (can't numerically sort symbols)
        boundaries = sorted(set(boundaries), key=str)

        # Build time regions:
        # Region 0: 0 < t < boundaries[0]
        # Region 1: boundaries[0] < t < boundaries[1]
        # ...
        # Region N: t > boundaries[-1]
        region_bounds: list[tuple[sp.Expr, sp.Expr]] = []
        prev = sp.S.Zero
        for bp in boundaries:
            region_bounds.append((prev, bp))
            prev = bp
        region_bounds.append((prev, oo))

        steps.append((
            "Piecewise forcing definition",
            "f(t) = {\n" + "\n".join(desc_lines) + "\n}\n\n"
            f"Boundary points: {', '.join(str(b) for b in boundaries)}\n"
            f"Time regions: {len(region_bounds)}",
        ))

        steps.append((
            "Transition matrix approach",
            "각 구간에서:\n"
            "  x(t) = Φ(t - tₖ₋₁)·x(tₖ₋₁) + ∫_{tₖ₋₁}^{t} Φ(t-τ)·B·fₖ(τ) dτ\n"
            "  f(t)=0인 구간: x(t) = Φ(t - tₖ₋₁)·x(tₖ₋₁)  (자유진동)",
        ))

        # ── Solve region by region ──
        x_current = x0_sym  # state at start of current region
        all_regions = []

        for region_idx, (r_lo, r_hi) in enumerate(region_bounds):
            region_num = region_idx + 1
            if r_hi == oo:
                region_label = f"Region {region_num} (t > {r_lo})"
            else:
                region_label = f"Region {region_num} ({r_lo} < t < {r_hi})"

            # Which piece is active in this region?
            f_active = sp.S.Zero
            for expr, p_lo, p_hi in parsed_pieces:
                # Check if region overlaps with piece interval
                if self._intervals_overlap(r_lo, r_hi, p_lo, p_hi):
                    f_active = expr
                    break

            if f_active == 0 or f_active == sp.S.Zero:
                # Pure free vibration: x(t) = Φ(t - r_lo) · x(r_lo)
                Phi_shift = Phi.subs(t_s, t_s - r_lo)
                x_region = simplify(Phi_shift * x_current)

                r_lines = [f"  x_{i + 1}(t) = {x_region[i, 0]}" for i in range(n)]
                steps.append((
                    region_label,
                    f"f(t) = 0  →  자유진동\n"
                    f"x(t) = Φ(t - {r_lo}) · x({r_lo})\n\n"
                    + "\n".join(r_lines),
                ))
            else:
                # Forced: x(t) = Φ(t-r_lo)·x(r_lo) + ∫_{r_lo}^{t} Φ(t-τ)·B·f(τ) dτ
                # Free part
                Phi_shift_free = Phi.subs(t_s, t_s - r_lo)
                x_free_part = Phi_shift_free * x_current

                # Forced part: ∫_{r_lo}^{t} Φ(t-τ)·B·f(τ) dτ
                integrand = PhiB * f_active
                x_forced_part = sp.zeros(n, 1)
                for i in range(n):
                    result = self._sym_integrate_conv(
                        integrand[i, 0], tau_s, r_lo, t_s,
                    )
                    x_forced_part[i, 0] = result

                x_region = sp.zeros(n, 1)
                r_lines = []
                for i in range(n):
                    total = simplify(x_free_part[i, 0] + x_forced_part[i, 0])
                    x_region[i, 0] = total
                    r_lines.append(f"  x_{i + 1}(t) = {total}")

                steps.append((
                    region_label,
                    f"f(τ) = {f_active}\n"
                    f"x(t) = Φ(t-{r_lo})·x({r_lo}) + ∫_{{{r_lo}}}^{{t}} Φ(t-τ)·B·f(τ) dτ\n\n"
                    + "\n".join(r_lines),
                ))

            # Store region result
            all_regions.append({
                "label": region_label,
                "x_expr": [str(x_region[i, 0]) for i in range(n)],
            })
            final_answer[f"x_region{region_num}"] = [
                str(x_region[i, 0]) for i in range(n)
            ]

            # Evaluate state at end of region (= start of next)
            if r_hi != oo and not r_hi.is_infinite:
                x_at_boundary = sp.zeros(n, 1)
                bnd_lines = []
                for i in range(n):
                    val = simplify(x_region[i, 0].subs(t_s, r_hi))
                    x_at_boundary[i, 0] = val
                    bnd_lines.append(f"  x_{i + 1}({r_hi}) = {val}")
                steps.append((
                    f"State at t = {r_hi}",
                    "\n".join(bnd_lines),
                ))
                final_answer[f"x_at_{r_hi}"] = [
                    str(x_at_boundary[i, 0]) for i in range(n)
                ]
                x_current = x_at_boundary

        # Final summary
        summary_lines = []
        for reg in all_regions:
            summary_lines.append(f"{reg['label']}:")
            for i, expr in enumerate(reg["x_expr"]):
                summary_lines.append(f"  x_{i + 1}(t) = {expr}")
        final_answer["forced_response_symbolic"] = "\n".join(summary_lines)

    @staticmethod
    def _intervals_overlap(a_lo, a_hi, b_lo, b_hi) -> bool:
        """Check if two symbolic intervals could overlap (conservative)."""
        import sympy as sp
        # For symbolic bounds, assume they overlap if ordering is consistent
        # Simple heuristic: same start or str-based ordering
        if a_lo == b_lo:
            return True
        if a_hi == sp.oo or b_hi == sp.oo:
            if str(a_lo) >= str(b_lo) and (b_hi == sp.oo or str(a_lo) < str(b_hi)):
                return True
            if str(b_lo) >= str(a_lo) and (a_hi == sp.oo or str(b_lo) < str(a_hi)):
                return True
        # Numeric comparison if possible
        try:
            a_lo_f, a_hi_f = float(a_lo), float(a_hi)
            b_lo_f, b_hi_f = float(b_lo), float(b_hi)
            return a_lo_f < b_hi_f and b_lo_f < a_hi_f
        except (TypeError, ValueError):
            pass
        # Fallback: check if starts match region
        return str(a_lo) == str(b_lo)

    # -------------------------------------------------------------------
    # Core analysis engine (numeric)
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

        if f_type in ("exponential", "exp_truncated"):
            self._symbolic_transition_convolution(
                A, B, x0, n, forcing, steps, final_answer,
            )
            return

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
    # Symbolic forced response via transition matrix convolution
    # -------------------------------------------------------------------
    def _symbolic_transition_convolution(
        self,
        A: np.ndarray,
        B: np.ndarray,
        x0: np.ndarray,
        n: int,
        forcing: dict,
        steps: list[tuple[str, Any]],
        final_answer: dict[str, Any],
    ) -> None:
        """전이행렬 기반 심볼릭 강제응답: x(t) = Φ(t)x₀ + ∫₀ᵗ Φ(t-τ)Bf(τ)dτ."""
        try:
            import sympy as sp
            from sympy import (
                Symbol, Matrix, eye, exp, cos, Rational,
                simplify, trigsimp, nsimplify,
                inverse_laplace_transform, Heaviside,
            )
        except ImportError:
            return

        t_s = Symbol("t", positive=True)
        tau_s = Symbol("tau", positive=True)
        s = Symbol("s")

        # ── symbolic matrices ──
        def _to_sym(val):
            return Rational(val).limit_denominator(10000)

        A_s = Matrix([[_to_sym(x) for x in row] for row in A.tolist()])
        B_s = Matrix([[_to_sym(x) for x in row] for row in B.tolist()])
        x0_s = Matrix([_to_sym(x) for x in x0.tolist()])

        # ── Φ(t) via Laplace: L⁻¹{(sI - A)⁻¹} ──
        sI_A = s * eye(n) - A_s
        char_poly = sp.expand(sI_A.det())
        adj = sI_A.adjugate()

        Phi = sp.zeros(n, n)
        for i in range(n):
            for j in range(n):
                F_s = adj[i, j] / char_poly
                f_t = inverse_laplace_transform(F_s, s, t_s)
                f_t = f_t.rewrite(sp.exp)
                f_t = f_t.replace(Heaviside(t_s), sp.Integer(1))
                f_t = f_t.replace(Heaviside, lambda *a: sp.Integer(1))
                f_t = simplify(trigsimp(f_t.rewrite(cos)))
                Phi[i, j] = f_t

        steps.append((
            "Transition matrix Φ(t) (for convolution)",
            "\n".join(
                f"  Φ_{i + 1}{j + 1}(t) = {Phi[i, j]}"
                for i in range(n) for j in range(n)
            ),
        ))

        # ── parse forcing parameters ──
        def _pp(val, name, **kw):
            if isinstance(val, str):
                return Symbol(name, **kw)
            return nsimplify(val)

        f_type = forcing.get("type", "")
        F0 = _pp(forcing.get("F0", 1), "F_0", positive=True)
        a_p = _pp(forcing.get("a", 1), "a", positive=True)

        # ── Φ(t-τ)·B ──
        Phi_shifted = Phi.subs(t_s, t_s - tau_s)
        PhiB = Phi_shifted * B_s  # (n × m), typically n × 1

        if f_type == "exponential":
            f_tau = F0 * exp(-a_p * tau_s)
            integrand = PhiB * f_tau

            steps.append((
                "Convolution integral (exponential forcing)",
                f"f(t) = {F0}·e^(-{a_p}·t)\n\n"
                f"x(t) = Φ(t)·x₀ + ∫₀ᵗ Φ(t-τ)·B·{F0}·e^(-{a_p}·τ) dτ",
            ))

            x_forced = sp.zeros(n, 1)
            r_lines = []
            for i in range(n):
                result = self._sym_integrate_conv(
                    integrand[i, 0], tau_s, sp.S.Zero, t_s,
                )
                x_forced[i, 0] = result

            x_free = Phi * x0_s
            x_total = sp.zeros(n, 1)
            for i in range(n):
                x_total[i, 0] = simplify(x_free[i, 0] + x_forced[i, 0])
                r_lines.append(f"  x_{i + 1}(t) = {x_total[i, 0]}")

            steps.append((
                "Closed-form response (exponential forcing)",
                "\n".join(r_lines),
            ))
            final_answer["forced_response_symbolic"] = "\n".join(r_lines)

        elif f_type == "exp_truncated":
            t1 = _pp(forcing.get("t1", 1), "t_1", positive=True)
            f_tau = F0 * exp(-a_p * tau_s)
            integrand = PhiB * f_tau

            steps.append((
                "Forced response via transition matrix (truncated exponential)",
                f"f(t) = {F0}·e^(-{a_p}·t),  0 < t < {t1}\n"
                f"f(t) = 0,                    t > {t1}\n\n"
                f"Zero ICs → x(t) = ∫₀ᵗ Φ(t-τ)·B·f(τ) dτ\n\n"
                f"Region I  (0 < t ≤ {t1}): upper limit = t\n"
                f"Region II (t > {t1}):      x(t) = Φ(t-{t1})·x({t1})  (free vibration)",
            ))

            # ── Region I: 0 < t ≤ t1 ──
            x_r1 = sp.zeros(n, 1)
            r1_lines = []
            for i in range(n):
                result = self._sym_integrate_conv(
                    integrand[i, 0], tau_s, sp.S.Zero, t_s,
                )
                x_r1[i, 0] = result
                r1_lines.append(f"  x_{i + 1}(t) = {result}")

            steps.append((f"Region I (0 < t ≤ {t1})", "\n".join(r1_lines)))

            # ── state at t = t1 ──
            x_t1 = sp.zeros(n, 1)
            t1_lines = []
            for i in range(n):
                val = simplify(x_r1[i, 0].subs(t_s, t1))
                x_t1[i, 0] = val
                t1_lines.append(f"  x_{i + 1}({t1}) = {val}")

            steps.append((f"State vector at t = {t1}", "\n".join(t1_lines)))

            # ── Region II: t > t1 ──
            Phi_tt1 = Phi.subs(t_s, t_s - t1)
            x_r2 = sp.zeros(n, 1)
            x_r2_raw = Phi_tt1 * x_t1
            r2_lines = []
            for i in range(n):
                val = simplify(x_r2_raw[i, 0])
                x_r2[i, 0] = val
                r2_lines.append(f"  x_{i + 1}(t) = {val}")

            steps.append((
                f"Region II (t > {t1})",
                f"x(t) = Φ(t - {t1}) · x({t1})\n\n" + "\n".join(r2_lines),
            ))

            final_answer["x_region1"] = [str(x_r1[i, 0]) for i in range(n)]
            final_answer["x_region2"] = [str(x_r2[i, 0]) for i in range(n)]
            final_answer["x_at_t1"] = [str(x_t1[i, 0]) for i in range(n)]
            final_answer["forced_response_symbolic"] = (
                f"Region I (0 < t ≤ {t1}):\n"
                + "\n".join(r1_lines)
                + f"\n\nRegion II (t > {t1}):\n"
                + "\n".join(r2_lines)
            )

    @staticmethod
    def _sym_integrate_conv(integrand, var, lo, hi):
        """Integrate by expanding trig to separate t and τ dependencies."""
        import sympy as sp

        expanded = sp.expand(sp.expand_trig(sp.expand(integrand)))
        terms = sp.Add.make_args(expanded)

        result = sp.S.Zero
        for term in terms:
            factors = sp.Mul.make_args(term)
            outer = sp.S.One
            inner = sp.S.One
            for f in factors:
                if f.has(var):
                    inner *= f
                else:
                    outer *= f
            int_result = sp.integrate(inner, (var, lo, hi))
            if isinstance(int_result, sp.Piecewise) and len(int_result.args) >= 1:
                int_result = int_result.args[0][0]
            result += outer * int_result

        return sp.simplify(result)

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
