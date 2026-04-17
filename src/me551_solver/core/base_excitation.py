"""Base Excitation (지반 가진) 솔버.

SDOF: Mx'' + c(x'-y') + k(x-y) = 0,  y(t) = A sin(wt)
MDOF: Mx'' + Cx' + Kx = Ci y' + Ki y,  y(t) = A sin(wt)

Produces step-by-step analytical derivation of magnitude and phase angle,
matching standard textbook/exam solution format.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np

from .base import BaseSolver, SolverResult


# ---------------------------------------------------------------------------
# Equation parser
# ---------------------------------------------------------------------------

def parse_base_excitation_eq(eom_str: str, excitation_str: str) -> dict:
    """Parse base excitation equation strings into parameter dict.

    Supports various input forms:

    EOM forms:
      "2x'' + 3x' + 5x = 3y' + 5y"
      "Mx'' + cx' + kx = cy' + ky"           (symbolic)
      "x'' + 2x' + 4x = 2y' + 4y"
      "5x'' + 10(x'-y') + 20(x-y) = 0"      (Newton form)
      "Mx'' + c(x'-y') + k(x-y) = 0"        (symbolic Newton)

    Excitation forms:
      "y = 0.05 sin 8t"
      "y = A sin wt"
      "y(t) = 3sin(10t)"
      "y = 2 sin(5t)"

    Returns dict with keys: m, c, k, A, omega (float or None for symbolic).
    """
    eom = eom_str.strip()
    exc = excitation_str.strip()

    result: dict[str, Any] = {}

    # --- Parse excitation: y = A sin(wt) or y = A cos(wt) ---
    exc_clean = exc.replace(" ", "")
    # Remove "y(t)=" or "y="
    exc_clean = re.sub(r"^y\(?t?\)?\s*=\s*", "", exc_clean, flags=re.IGNORECASE)

    # Match: [coeff] (sin|cos) [( ] [freq] [*] t [) ]
    # Strategy: find sin/cos first, then split into amplitude (before) and freq*t (after)
    m_trig = re.search(r"(sin|cos)", exc_clean, re.IGNORECASE)
    if m_trig:
        func_str = m_trig.group(1).lower()
        amp_str = exc_clean[:m_trig.start()].strip()
        after = exc_clean[m_trig.end():].strip()

        # Parse amplitude
        if amp_str in ("", "+"):
            result["A"] = 1.0
        elif amp_str == "-":
            result["A"] = -1.0
        else:
            try:
                result["A"] = float(amp_str)
            except ValueError:
                result["A"] = None  # symbolic (e.g., "A")

        # Parse freq*t from after part: "8t", "(8t)", "wt", "(w*t)", "(omega*t)"
        after = after.strip("() ")
        # Remove trailing "t" and optional "*" before it
        m_freq = re.match(r"(.*?)\s*\*?\s*t\s*$", after, re.IGNORECASE)
        if m_freq:
            freq_str = m_freq.group(1).strip()
            if freq_str == "":
                result["omega"] = 1.0  # just "t" means w=1
            else:
                try:
                    result["omega"] = float(freq_str)
                except ValueError:
                    result["omega"] = None  # symbolic (e.g., "w", "omega")
        else:
            result["omega"] = None

        result["func"] = func_str
    else:
        result["A"] = None
        result["omega"] = None
        result["func"] = "sin"  # default

    # --- Parse EOM ---
    eom_clean = eom.replace(" ", "").replace("−", "-").replace("–", "-")

    # --- Try Newton form first: ax'' + b(x'-y') + c(x-y) = 0 ---
    m_newton = re.match(
        r"([a-zA-Z_]\w*|\d*\.?\d*)\s*x['\u2032]{2}\s*"    # M x''
        r"\+\s*([a-zA-Z_]\w*|\d*\.?\d*)\s*"                # + c
        r"\(\s*x['\u2032]\s*-\s*y['\u2032]\s*\)\s*"        # (x'-y')
        r"\+\s*([a-zA-Z_]\w*|\d*\.?\d*)\s*"                # + k
        r"\(\s*x\s*-\s*y\s*\)\s*=\s*0$",                   # (x-y)=0
        eom_clean, re.IGNORECASE,
    )
    if m_newton:
        for i, key in enumerate(["m", "c", "k"]):
            val_str = m_newton.group(i + 1)
            if not val_str or val_str == "":
                result[key] = 1.0
            else:
                try:
                    result[key] = float(val_str)
                except ValueError:
                    result[key] = None  # symbolic
        return result

    # --- Standard form: ax'' + bx' + cx = dy' + ey ---
    # Split on "="
    if "=" in eom_clean:
        lhs, rhs = eom_clean.split("=", 1)
    else:
        lhs = eom_clean
        rhs = ""

    _COEFF_PAT = r"([+-]?(?:\d+\.?\d*|[a-zA-Z_]\w*)?)"

    _NOT_FOUND = object()  # sentinel

    def _extract_coeffs(expr: str, var: str) -> dict[str, Any]:
        """Extract coefficients of var'', var', var from expression.
        Returns _NOT_FOUND if term not found, None if symbolic, float if numeric.
        """
        coeffs: dict[str, Any] = {"dd": _NOT_FOUND, "d": _NOT_FOUND, "0": _NOT_FOUND}
        CP = _COEFF_PAT

        def _parse_coeff(cs: str) -> float | None:
            cs = cs.strip()
            if cs in ("", "+"):
                return 1.0
            if cs == "-":
                return -1.0
            # Strip leading + for things like "+c" or "+3"
            if cs.startswith("+"):
                cs = cs[1:]
            try:
                return float(cs)
            except ValueError:
                return None  # symbolic like "M", "c", "k"

        # Double derivative: coeff * var''
        for pat_dd in [
            CP + rf"\s*\*?\s*{var}['\u2032]{{2}}",
            CP + rf"\s*\*?\s*{var}''",
        ]:
            m_dd = re.search(pat_dd, expr)
            if m_dd:
                coeffs["dd"] = _parse_coeff(m_dd.group(1))
                break

        # Single derivative: coeff * var'  (but not var'')
        # Remove double-deriv matches first
        expr_no_dd = re.sub(r"[a-zA-Z_\w\d.]*\s*\*?\s*" + var + r"['\u2032]{2}", "", expr)
        expr_no_dd = re.sub(r"[a-zA-Z_\w\d.]*\s*\*?\s*" + var + r"''", "", expr_no_dd)

        for pat_d in [
            CP + rf"\s*\*?\s*{var}['\u2032]",
            CP + rf"\s*\*?\s*{var}'",
        ]:
            m_d = re.search(pat_d, expr_no_dd)
            if m_d:
                coeffs["d"] = _parse_coeff(m_d.group(1))
                break

        # Zeroth order: coeff * var  (but not var' or var'')
        expr_no_deriv = re.sub(r"[a-zA-Z_\w\d.]*\s*\*?\s*" + var + r"['\u2032]{1,2}", "", expr_no_dd)
        expr_no_deriv = re.sub(r"[a-zA-Z_\w\d.]*\s*\*?\s*" + var + r"'{1,2}", "", expr_no_deriv)

        m_0 = re.search(CP + rf"\s*\*?\s*{var}(?!['\u2032a-zA-Z])", expr_no_deriv)
        if m_0:
            coeffs["0"] = _parse_coeff(m_0.group(1))

        return coeffs

    x_coeffs = _extract_coeffs(lhs, "x")
    y_coeffs = _extract_coeffs(rhs, "y")

    # M = coeff of x''
    if x_coeffs["dd"] is not _NOT_FOUND:
        result["m"] = x_coeffs["dd"]  # float or None (symbolic)

    # c = coeff of x' (LHS) or y' (RHS)
    if x_coeffs["d"] is not _NOT_FOUND:
        result["c"] = x_coeffs["d"]
    elif y_coeffs["d"] is not _NOT_FOUND:
        result["c"] = y_coeffs["d"]

    # k = coeff of x (LHS) or y (RHS)
    if x_coeffs["0"] is not _NOT_FOUND:
        result["k"] = x_coeffs["0"]
    elif y_coeffs["0"] is not _NOT_FOUND:
        result["k"] = y_coeffs["0"]

    return result


class BaseExcitationSolver(BaseSolver):
    """Base excitation 해석 솔버 -- 수식 유도 + 수치 평가."""

    def solve(self, params: dict) -> SolverResult:
        n_dof = params.get("n_dof", 1)
        if n_dof == 1:
            return self._solve_sdof(params)
        else:
            return self._solve_mdof(params)

    # ===================================================================
    # SDOF -- step-by-step derivation
    # ===================================================================
    def _solve_sdof(self, params: dict) -> SolverResult:
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # --- Optional numeric values (None = symbolic) ---
        m = params.get("m")
        c = params.get("c")
        k = params.get("k")
        A = params.get("A")
        omega = params.get("omega")
        func = params.get("func", "sin")  # "sin" or "cos"
        has_numeric = all(v is not None for v in [m, c, k, A, omega])

        is_cos = (func == "cos")

        # ---------------------------------------------------------------
        # Step 1: Equation of motion
        # ---------------------------------------------------------------
        if is_cos:
            steps.append((
                "Step 1: Equation of Motion",
                "The support undergoes harmonic motion y(t) = A cos(wt).\n\n"
                "From Newton's 2nd law (or free-body diagram):\n"
                "  Mx'' + c(x' - y') + k(x - y) = 0\n\n"
                "Rearranging:\n"
                "  Mx'' + cx' + kx = cy' + ky\n\n"
                "Substituting y(t) = A cos(wt),  y'(t) = -Aw sin(wt):\n"
                "  Mx'' + cx' + kx = -Acw sin(wt) + Ak cos(wt)"
            ))
        else:
            steps.append((
                "Step 1: Equation of Motion",
                "The support undergoes harmonic motion y(t) = A sin(wt).\n\n"
                "From Newton's 2nd law (or free-body diagram):\n"
                "  Mx'' + c(x' - y') + k(x - y) = 0\n\n"
                "Rearranging:\n"
                "  Mx'' + cx' + kx = cy' + ky\n\n"
                "Substituting y(t) = A sin(wt),  y'(t) = Aw cos(wt):\n"
                "  Mx'' + cx' + kx = Acw cos(wt) + Ak sin(wt)"
            ))

        # ---------------------------------------------------------------
        # Step 2: Normalize by M
        # ---------------------------------------------------------------
        if is_cos:
            steps.append((
                "Step 2: Normalize (divide by M)",
                "Letting  w_n = sqrt(k/M),   zeta = c / (2*M*w_n),   r = w/w_n :\n\n"
                "  x'' + 2*zeta*w_n*x' + w_n^2*x = -2*zeta*w_n*A*w*sin(wt) + A*w_n^2*cos(wt)"
            ))
        else:
            steps.append((
                "Step 2: Normalize (divide by M)",
                "Letting  w_n = sqrt(k/M),   zeta = c / (2*M*w_n),   r = w/w_n :\n\n"
                "  x'' + 2*zeta*w_n*x' + w_n^2*x = 2*zeta*w_n*A*w*cos(wt) + A*w_n^2*sin(wt)"
            ))

        # ---------------------------------------------------------------
        # Step 3: Combine RHS into single harmonic
        # ---------------------------------------------------------------
        if is_cos:
            steps.append((
                "Step 3: Combine RHS into single harmonic",
                "The RHS has the form  P*cos(wt) + Q*sin(wt)  where:\n"
                "  P = A*w_n^2\n"
                "  Q = -2*zeta*w_n*A*w = -A*w_n^2 * 2*zeta*r\n\n"
                "Combined amplitude:\n"
                "  sqrt(P^2 + Q^2) = A*w_n^2 * sqrt(1 + (2*zeta*r)^2)\n\n"
                "Therefore:\n"
                "  RHS = A*w_n^2 * sqrt(1 + (2*zeta*r)^2) * cos(wt + alpha)\n\n"
                "where:\n"
                "  tan(alpha) = -Q/P = 2*zeta*r"
            ))
        else:
            steps.append((
                "Step 3: Combine RHS into single harmonic",
                "The RHS has the form  P*cos(wt) + Q*sin(wt)  where:\n"
                "  P = 2*zeta*w_n*A*w = A*w_n^2 * 2*zeta*r\n"
                "  Q = A*w_n^2\n\n"
                "Combined amplitude:\n"
                "  sqrt(P^2 + Q^2) = A*w_n^2 * sqrt((2*zeta*r)^2 + 1)\n"
                "                  = A*w_n^2 * sqrt(1 + (2*zeta*r)^2)\n\n"
                "Therefore:\n"
                "  RHS = A*w_n^2 * sqrt(1 + (2*zeta*r)^2) * cos(wt - alpha)\n\n"
                "where:\n"
                "  tan(alpha) = Q/P = (A*w_n^2) / (A*w_n^2 * 2*zeta*r) = 1 / (2*zeta*r)"
            ))

        # ---------------------------------------------------------------
        # Step 4: Apply frequency response function G(jw)
        # ---------------------------------------------------------------
        _alpha_sign = "+" if is_cos else "-"
        steps.append((
            "Step 4: System response via G(jw)",
            "The normalized EOM is:\n"
            f"  x'' + 2*zeta*w_n*x' + w_n^2*x = F_0 * cos(wt {_alpha_sign} alpha)\n\n"
            "where F_0 = A*w_n^2 * sqrt(1 + (2*zeta*r)^2).\n\n"
            "The frequency response function of the standard 2nd-order system is:\n"
            "  G(jw) = 1 / (w_n^2 - w^2 + j*2*zeta*w_n*w)\n"
            "        = (1/w_n^2) * 1 / (1 - r^2 + j*2*zeta*r)\n\n"
            "|G(jw)| = (1/w_n^2) / sqrt((1 - r^2)^2 + (2*zeta*r)^2)\n\n"
            "The particular solution has the form:\n"
            f"  x(t) = F_0 * |G(jw)| * cos(wt {_alpha_sign} alpha - phi)"
        ))

        # ---------------------------------------------------------------
        # Step 5: Phase angle phi from G(jw)
        # ---------------------------------------------------------------
        steps.append((
            "Step 5: Phase angle phi of G(jw)",
            "From G(jw) = 1 / (1 - r^2 + j*2*zeta*r):\n\n"
            "  tan(phi) = -Im(G) / Re(G) = 2*zeta*r / (1 - r^2)\n\n"
            "  phi = arctan(2*zeta*r / (1 - r^2))"
        ))

        # ---------------------------------------------------------------
        # Step 6: Final magnitude (same for sin and cos)
        # ---------------------------------------------------------------
        steps.append((
            "Step 6: Magnitude of steady-state response",
            "Magnitude = F_0 * |G(jw)|\n"
            "          = A*w_n^2 * sqrt(1 + (2*zeta*r)^2) * (1/w_n^2) / sqrt((1-r^2)^2 + (2*zeta*r)^2)\n\n"
            "  +-------------------------------------------------+\n"
            "  |                                                 |\n"
            "  |              A * sqrt(1 + (2*zeta*r)^2)         |\n"
            "  |  |X| =  ---------------------------------      |\n"
            "  |          sqrt((1 - r^2)^2 + (2*zeta*r)^2)      |\n"
            "  |                                                 |\n"
            "  +-------------------------------------------------+"
        ))

        # ---------------------------------------------------------------
        # Step 7: Total phase angle
        # ---------------------------------------------------------------
        if is_cos:
            steps.append((
                "Step 7: Total phase angle",
                "Total phase = -alpha + phi  (note: +alpha in forcing becomes -alpha after G(jw))\n\n"
                "  alpha = arctan(2*zeta*r)\n"
                "  phi   = arctan(2*zeta*r / (1 - r^2))\n\n"
                "  +------------------------------------------------------+\n"
                "  |                                                      |\n"
                "  |  Total phase = phi - alpha                           |\n"
                "  |    = arctan(2*zeta*r / (1-r^2)) - arctan(2*zeta*r)  |\n"
                "  |                                                      |\n"
                "  +------------------------------------------------------+"
            ))
        else:
            steps.append((
                "Step 7: Total phase angle",
                "Total phase = alpha + phi\n\n"
                "  alpha = arctan(1 / (2*zeta*r))\n"
                "  phi   = arctan(2*zeta*r / (1 - r^2))\n\n"
                "Combining using tan(alpha + phi):\n\n"
                "  tan(alpha + phi) = (tan(alpha) + tan(phi)) / (1 - tan(alpha)*tan(phi))\n\n"
                "  Numerator:  1/(2*zeta*r) + 2*zeta*r/(1-r^2)\n"
                "            = [(1-r^2) + (2*zeta*r)^2] / [2*zeta*r * (1-r^2)]\n\n"
                "  Denominator: 1 - [1/(2*zeta*r)] * [2*zeta*r/(1-r^2)]\n"
                "             = 1 - 1/(1-r^2)\n"
                "             = -r^2 / (1-r^2)\n\n"
                "  tan(alpha+phi) = [(1-r^2) + (2*zeta*r)^2] / [2*zeta*r*(1-r^2)]\n"
                "                   * [(1-r^2) / (-r^2)]\n\n"
                "  +----------------------------------------------------------+\n"
                "  |                                                          |\n"
                "  |                  -[ (2*zeta*r)^2 + 1 - r^2 ]            |\n"
                "  |  alpha + phi = arctan --------------------------        |\n"
                "  |                          2 * zeta * r^3                  |\n"
                "  |                                                          |\n"
                "  +----------------------------------------------------------+"
            ))

        # ---------------------------------------------------------------
        # Step 8: Final answer (symbolic)
        # ---------------------------------------------------------------
        if is_cos:
            steps.append((
                "Step 8: Complete solution",
                "x(t) = |X| * cos(wt - (phi - alpha))\n\n"
                "where:\n\n"
                "  (Magnitude)    |X| = A * sqrt(1 + (2*zeta*r)^2) / sqrt((1-r^2)^2 + (2*zeta*r)^2)\n\n"
                "  (Phase angle)  phi - alpha = arctan(2*zeta*r / (1-r^2)) - arctan(2*zeta*r)\n\n"
                "  with  r = w/w_n,   zeta = c/(2*M*w_n),   w_n = sqrt(k/M)"
            ))
            final_answer["magnitude_formula"] = (
                "|X| = A * sqrt(1 + (2*zeta*r)^2) / sqrt((1-r^2)^2 + (2*zeta*r)^2)"
            )
            final_answer["phase_formula"] = (
                "phi - alpha = arctan(2*zeta*r / (1-r^2)) - arctan(2*zeta*r)"
            )
        else:
            steps.append((
                "Step 8: Complete solution",
                "x(t) = |X| * cos(wt - alpha - phi)\n\n"
                "where:\n\n"
                "  (Magnitude)    |X| = A * sqrt(1 + (2*zeta*r)^2) / sqrt((1-r^2)^2 + (2*zeta*r)^2)\n\n"
                "  (Phase angle)  alpha + phi = arctan[ -((2*zeta*r)^2 + 1 - r^2) / (2*zeta*r^3) ]\n\n"
                "  with  r = w/w_n,   zeta = c/(2*M*w_n),   w_n = sqrt(k/M)"
            ))
            final_answer["magnitude_formula"] = (
                "|X| = A * sqrt(1 + (2*zeta*r)^2) / sqrt((1-r^2)^2 + (2*zeta*r)^2)"
            )
            final_answer["phase_formula"] = (
                "alpha + phi = arctan[ -((2*zeta*r)^2 + 1 - r^2) / (2*zeta*r^3) ]"
            )

        # ---------------------------------------------------------------
        # Step 9: Numeric evaluation (if values given)
        # ---------------------------------------------------------------
        if has_numeric:
            m = float(m)
            c = float(c)
            k = float(k)
            A = float(A)
            omega = float(omega)

            omega_n = np.sqrt(k / m)
            zeta = c / (2.0 * m * omega_n)
            r = omega / omega_n

            mag_num = 1.0 + (2.0 * zeta * r) ** 2
            mag_den = (1.0 - r**2) ** 2 + (2.0 * zeta * r) ** 2
            X_amp = A * np.sqrt(mag_num / mag_den)

            if is_cos:
                alpha = np.arctan(2.0 * zeta * r)
                phi = np.arctan2(2.0 * zeta * r, 1.0 - r**2)
                total_phase = phi - alpha
                alpha_formula = f"alpha = arctan({2*zeta*r:.6g}) = {np.degrees(alpha):.4f} deg"
                total_label = "phi - alpha"
            else:
                alpha = np.arctan2(1.0, 2.0 * zeta * r)
                phi = np.arctan2(2.0 * zeta * r, 1.0 - r**2)
                total_phase = alpha + phi
                alpha_formula = f"alpha = arctan(1 / {2*zeta*r:.6g}) = {np.degrees(alpha):.4f} deg"
                total_label = "alpha + phi"

            steps.append((
                "Step 9: Numeric evaluation",
                f"Given:  M = {m},  c = {c},  k = {k},  A = {A},  w = {omega}\n"
                f"        y(t) = {A} {func}({omega}t)\n\n"
                f"  w_n  = sqrt({k}/{m}) = {omega_n:.6g} rad/s\n"
                f"  zeta = {c}/(2*{m}*{omega_n:.6g}) = {zeta:.6g}\n"
                f"  r    = {omega}/{omega_n:.6g} = {r:.6g}\n"
                f"  2*zeta*r = {2*zeta*r:.6g}\n\n"
                f"Magnitude:\n"
                f"  |X| = {A} * sqrt(1 + {(2*zeta*r)**2:.6g}) / sqrt({(1-r**2)**2:.6g} + {(2*zeta*r)**2:.6g})\n"
                f"      = {A} * {np.sqrt(mag_num):.6g} / {np.sqrt(mag_den):.6g}\n"
                f"      = {X_amp:.6g}\n\n"
                f"Phase angle:\n"
                f"  {alpha_formula}\n"
                f"  phi   = arctan({2*zeta*r:.6g} / {1-r**2:.6g}) = {np.degrees(phi):.4f} deg\n"
                f"  {total_label} = {np.degrees(total_phase):.4f} deg\n\n"
                f"Response:\n"
                f"  x(t) = {X_amp:.6g} cos({omega}t - {total_phase:.6g})"
            ))

            final_answer["omega_n"] = omega_n
            final_answer["zeta"] = zeta
            final_answer["r"] = r
            final_answer["X_amplitude"] = X_amp
            final_answer["alpha_deg"] = np.degrees(alpha)
            final_answer["phi_deg"] = np.degrees(phi)
            final_answer["total_phase_deg"] = np.degrees(total_phase)

            given = {"M": m, "c": c, "k": k, "A": A, "omega": omega, "func": func}
        else:
            given = {"M": "M", "c": "c", "k": "k", "A": "A", "omega": "w", "func": func}

        # ---------------------------------------------------------------
        # Sanity checks
        # ---------------------------------------------------------------
        sanity_parts = []
        sanity_parts.append(
            "At r -> 0:  |X| -> A  (system follows base motion quasi-statically)"
        )
        sanity_parts.append(
            "At r = sqrt(2):  All T_d curves pass through T_d = 1  (crossover point)"
        )
        sanity_parts.append(
            "At r >> 1:  |X| -> 0  (isolation region, mass stays still)"
        )
        if has_numeric:
            T_d = X_amp / A
            if T_d > 1:
                sanity_parts.append(f"T_d = {T_d:.4g} > 1 : amplification region (r < sqrt(2))")
            else:
                sanity_parts.append(f"T_d = {T_d:.4g} <= 1 : isolation region (r > sqrt(2))")

        return SolverResult(
            problem_type="Base Excitation (SDOF)",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    # ===================================================================
    # MDOF
    # ===================================================================
    def _solve_mdof(self, params: dict) -> SolverResult:
        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}
        sanity_parts: list[str] = []

        M = np.array(params["M"], dtype=float)
        K = np.array(params["K"], dtype=float)
        n = M.shape[0]

        C_input = params.get("C")
        C = np.array(C_input, dtype=float) if C_input is not None else np.zeros_like(M)

        iota_input = params.get("influence_vector")
        iota = np.array(iota_input, dtype=float) if iota_input is not None else np.ones(n)

        A_amp = float(params.get("A", 1.0))
        omega = float(params.get("omega", 1.0))

        def _fmt_mat(mat, name):
            lines = [f"{name} ="]
            for row in mat:
                lines.append("  [" + ", ".join(f"{v:.6g}" for v in row) + "]")
            return "\n".join(lines)

        # Step 1: EOM
        steps.append((
            "Step 1: Equation of Motion",
            "Mx'' + C(x' - i*y') + K(x - i*y) = 0\n"
            "  => Mx'' + Cx' + Kx = Ci*y' + Ki*y\n\n"
            f"y(t) = {A_amp} sin({omega} t),  DOF = {n}\n"
            f"Influence vector i = {iota.tolist()}\n\n"
            + "\n\n".join([_fmt_mat(M, "M"), _fmt_mat(C, "C"), _fmt_mat(K, "K")])
        ))

        # Step 2: Relative coordinate
        F_eff = M @ iota * A_amp * omega**2
        steps.append((
            "Step 2: Relative coordinate z = x - i*y",
            "Mz'' + Cz' + Kz = -Mi*y'' = Mi*A*w^2 * sin(wt)\n\n"
            f"Effective force: Mi*A*w^2 = {F_eff.tolist()}"
        ))
        final_answer["F_eff"] = F_eff.tolist()

        # Step 3: Frequency response
        H_matrix = -omega**2 * M + 1j * omega * C + K
        try:
            H_inv = np.linalg.inv(H_matrix)
            Z_complex = H_inv @ F_eff
            Z_amp = np.abs(Z_complex)
            Z_phase = np.angle(Z_complex)

            X_complex = Z_complex + iota * A_amp
            X_amp = np.abs(X_complex)
            X_phase = np.angle(X_complex)

            z_lines = [f"  |Z_{i+1}| = {Z_amp[i]:.6g},  phase = {np.degrees(Z_phase[i]):.4f} deg" for i in range(n)]
            x_lines = [f"  |X_{i+1}| = {X_amp[i]:.6g},  phase = {np.degrees(X_phase[i]):.4f} deg" for i in range(n)]
            T_d_vec = X_amp / A_amp
            td_lines = [f"  T_d_{i+1} = |X_{i+1}|/A = {T_d_vec[i]:.6g}" for i in range(n)]

            steps.append((
                "Step 3: Steady-state response via impedance matrix",
                "Z(w) = [-w^2*M + j*w*C + K]^(-1) * Mi*A*w^2\n"
                "X = Z + i*A\n\n"
                "Relative displacement:\n" + "\n".join(z_lines) + "\n\n"
                "Absolute displacement:\n" + "\n".join(x_lines) + "\n\n"
                "Transmissibility:\n" + "\n".join(td_lines)
            ))

            final_answer["Z_amplitudes"] = Z_amp.tolist()
            final_answer["Z_phases_deg"] = np.degrees(Z_phase).tolist()
            final_answer["X_amplitudes"] = X_amp.tolist()
            final_answer["X_phases_deg"] = np.degrees(X_phase).tolist()
            final_answer["T_d"] = T_d_vec.tolist()

        except np.linalg.LinAlgError:
            steps.append(("Step 3", "Impedance matrix singular at this frequency."))

        # Sanity
        sanity_parts.append(f"M symmetric: {'PASS' if np.allclose(M, M.T) else 'FAIL'}")
        sanity_parts.append(f"K symmetric: {'PASS' if np.allclose(K, K.T) else 'FAIL'}")
        eig_M = np.linalg.eigvalsh(M)
        sanity_parts.append(f"M positive definite: {'PASS' if np.all(eig_M > 0) else 'FAIL'}")

        given = {
            "M": M.tolist(), "C": C.tolist(), "K": K.tolist(),
            "influence_vector": iota.tolist(),
            "A": A_amp, "omega": omega, "n_dof": n,
        }

        return SolverResult(
            problem_type="Base Excitation (MDOF)",
            given=given, steps=steps, final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    # -------------------------------------------------------------------
    # Input template
    # -------------------------------------------------------------------
    def get_input_template(self) -> dict:
        return {
            "n_dof": {
                "type": "int", "required": False, "default": 1,
                "description": "자유도. 1이면 SDOF, 2 이상이면 MDOF.",
            },
            "m": {"type": "float", "description": "질량 (SDOF). 생략시 symbolic."},
            "c": {"type": "float", "description": "감쇠 계수 (SDOF). 생략시 symbolic."},
            "k": {"type": "float", "description": "강성 (SDOF). 생략시 symbolic."},
            "A": {"type": "float", "description": "Base 진폭. 생략시 symbolic."},
            "omega": {"type": "float", "description": "가진 각진동수."},
            "M": {"type": "list[list[float]]", "description": "질량행렬 (MDOF)."},
            "C": {"type": "list[list[float]]", "description": "감쇠행렬 (MDOF)."},
            "K": {"type": "list[list[float]]", "description": "강성행렬 (MDOF)."},
            "influence_vector": {"type": "list[float]", "description": "영향벡터 (MDOF)."},
        }
