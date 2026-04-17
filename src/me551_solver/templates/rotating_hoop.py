"""회전 후프(Rotating Hoop) 문제 템플릿.

Bead of mass m sliding on a circular hoop of radius R rotating at angular
velocity Omega about the vertical axis.

Physics:
    Generalised coordinate: theta (angle from the bottom of the hoop)

    Kinetic energy:
        T = (1/2) m (R^2 * thetadot^2 + R^2 * Omega^2 * sin^2(theta))

    Potential energy (taking bottom as reference, V = -mgR cos(theta)):
        V = -m*g*R*cos(theta)

    Lagrangian: L = T - V

    EOM (Lagrange's equation):
        m*R^2 * theta'' + m*g*R*sin(theta) - m*R^2*Omega^2*sin(theta)*cos(theta) = 0
        → theta'' + (g/R)*sin(theta) - Omega^2*sin(theta)*cos(theta) = 0

Equilibria:
    1. theta = 0  (bottom)
    2. theta = pi  (top, always unstable — not analysed by default)
    3. cos(theta_eq) = g / (R * Omega^2)  (non-trivial, exists when Omega^2 > g/R)

Stability:
    Linearise theta = theta_eq + epsilon:
        theta_eq = 0:
            EOM → epsilon'' + (g/R - Omega^2) * epsilon = 0
            STABLE  iff Omega < sqrt(g/R)
            UNSTABLE iff Omega > sqrt(g/R)

        theta_eq = arccos(g/(R*Omega^2)):
            EOM → epsilon'' + Omega^2 * sin^2(theta_eq) * epsilon = 0
            omega_n = Omega * |sin(theta_eq)|
            Always STABLE (for Omega > sqrt(g/R), coefficient > 0)
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from ..core.base import BaseSolver, SolverResult

# ---------------------------------------------------------------------------
# Module-level symbols (reused across methods)
# ---------------------------------------------------------------------------
_theta, _t = sp.symbols("theta t")
_m, _R, _Omega, _g = sp.symbols("m R Omega g", positive=True)
_eps = sp.symbols("epsilon")


class RotatingHoopSolver(BaseSolver):
    """회전 후프 위의 비드 운동: Lagrange → EOM → 평형점 → 선형화 → 안정성.

    Supports both fully symbolic (None parameters) and numeric inputs.
    """

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self, params: dict) -> SolverResult:
        """Analyse the rotating hoop system.

        Parameters
        ----------
        params : dict
            R        : float | None  – hoop radius  (symbolic if None)
            m        : float | None  – bead mass    (symbolic if None)
            Omega    : float | None  – angular velocity (symbolic if None)
            g        : float | None  – gravitational acceleration
                                       (symbolic if None; default symbol g)
            equilibrium_to_analyze : str
                "all"           – analyse both equilibria (default)
                "theta=0"       – analyse only the trivial equilibrium
                "theta=arccos"  – analyse only the non-trivial equilibrium
        """
        R_val     = params.get("R")
        m_val     = params.get("m")
        Omega_val = params.get("Omega")
        g_val     = params.get("g")
        eq_choice = params.get("equilibrium_to_analyze", "all")

        # Resolve symbolic / numeric parameters
        # Use sp.nsimplify for floats (avoids huge Rational representations)
        R_sym     = sp.nsimplify(R_val)     if R_val     is not None else _R
        m_sym     = sp.nsimplify(m_val)     if m_val     is not None else _m
        Omega_sym = sp.nsimplify(Omega_val) if Omega_val is not None else _Omega
        g_sym     = sp.nsimplify(g_val)     if g_val     is not None else _g

        given: dict[str, Any] = {
            "R (hoop radius)":        R_val  if R_val  is not None else "symbolic R",
            "m (bead mass)":          m_val  if m_val  is not None else "symbolic m",
            "Omega (angular vel.)":   Omega_val if Omega_val is not None else "symbolic Omega",
            "g (gravity)":            g_val  if g_val  is not None else "symbolic g",
            "Equilibria analysed":    eq_choice,
        }

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # Step 1 – Kinetic & potential energy
        theta = _theta
        T = sp.Rational(1, 2) * m_sym * (
            R_sym**2 * sp.Symbol("thetadot")**2
            + R_sym**2 * Omega_sym**2 * sp.sin(theta)**2
        )
        V = -m_sym * g_sym * R_sym * sp.cos(theta)
        steps.append((
            "Kinetic energy T",
            f"T = (1/2)*m*(R²·θ̇² + R²·Ω²·sin²θ)\n"
            f"  = {sp.sstr(T)}"
        ))
        steps.append((
            "Potential energy V",
            f"V = -m·g·R·cos(θ)\n"
            f"  = {sp.sstr(V)}"
        ))

        # Step 2 – Lagrangian and EOM
        # L = T - V
        # d/dt(∂L/∂θ̇) - ∂L/∂θ = 0
        # m*R^2*theta'' + m*g*R*sin(theta) - m*R^2*Omega^2*sin(theta)*cos(theta) = 0
        # Divide through by m*R^2:
        eom_lhs = (g_sym / R_sym) * sp.sin(theta) - Omega_sym**2 * sp.sin(theta) * sp.cos(theta)
        steps.append((
            "Lagrange's equation (EOM) — divide by m·R²",
            f"θ'' + (g/R)·sin(θ) − Ω²·sin(θ)·cos(θ) = 0\n"
            f"θ'' + sin(θ)·[(g/R) − Ω²·cos(θ)] = 0\n"
            f"Symbolic form:\n"
            f"  theta'' + {sp.sstr(eom_lhs)} = 0"
        ))
        final_answer["EOM"] = "theta'' + (g/R)*sin(theta) - Omega^2*sin(theta)*cos(theta) = 0"

        # Step 3 – Equilibria
        # sin(theta_eq) = 0  → theta_eq = 0  (or pi, top)
        # sin(theta_eq) ≠ 0  → cos(theta_eq) = g/(R*Omega^2)
        steps.append((
            "Equilibrium conditions: set θ'' = 0 → sin(θ)·[(g/R) − Ω²cos(θ)] = 0",
            "Case A: sin(θ_eq) = 0  → θ_eq = 0  (bottom)\n"
            "Case B: cos(θ_eq) = g/(R·Ω²)  — non-trivial (exists if |g/(R·Ω²)| ≤ 1)"
        ))

        cos_eq2 = g_sym / (R_sym * Omega_sym**2)
        cos_eq2_simplified = sp.simplify(cos_eq2)
        final_answer["equilibrium_theta0"] = "theta = 0 (always)"
        final_answer["equilibrium_nontrivial_cos"] = str(cos_eq2_simplified)

        # Numeric value of cos(theta_eq) if fully specified
        if all(v is not None for v in [R_val, Omega_val, g_val]):
            cos_eq2_num = float(cos_eq2_simplified)
            if abs(cos_eq2_num) <= 1.0:
                import math
                theta_eq2_num = math.acos(cos_eq2_num)
                final_answer["equilibrium_cos_theta"] = round(cos_eq2_num, 8)
                final_answer["equilibrium_theta_rad"] = round(theta_eq2_num, 8)
                steps.append((
                    "Non-trivial equilibrium (numeric)",
                    f"cos(θ_eq) = g/(R·Ω²) = {cos_eq2_num:.6f}\n"
                    f"θ_eq = arccos({cos_eq2_num:.6f}) = {theta_eq2_num:.6f} rad"
                    f" ({math.degrees(theta_eq2_num):.4f}°)"
                ))
            else:
                steps.append((
                    "Non-trivial equilibrium (numeric)",
                    f"cos(θ_eq) = g/(R·Ω²) = {cos_eq2_num:.6f} — outside [-1,1].\n"
                    f"Non-trivial equilibrium does NOT exist for these parameters."
                ))
                final_answer["equilibrium_nontrivial_exists"] = False
        else:
            steps.append((
                "Non-trivial equilibrium (symbolic)",
                f"cos(θ_eq) = g/(R·Ω²) = {sp.sstr(cos_eq2_simplified)}\n"
                f"Exists when Ω² > g/R  (i.e., Ω > √(g/R))"
            ))

        # Step 4 – Linearise about theta = 0
        if eq_choice in ("all", "theta=0"):
            steps.append((
                "Linearisation about θ_eq = 0",
                "Let θ = 0 + ε  (ε small)\n"
                "sin(ε) ≈ ε,  cos(ε) ≈ 1\n"
                "EOM → ε'' + [(g/R) − Ω²]·ε = 0"
            ))
            coeff_eq0 = g_sym / R_sym - Omega_sym**2
            coeff_eq0_simplified = sp.simplify(coeff_eq0)
            steps.append((
                "Linearised coefficient at θ=0",
                f"κ₀ = g/R − Ω²  =  {sp.sstr(coeff_eq0_simplified)}\n"
                f"  κ₀ > 0  →  STABLE   (Ω < √(g/R))\n"
                f"  κ₀ < 0  →  UNSTABLE (Ω > √(g/R))\n"
                f"  κ₀ = 0  →  MARGINALLY STABLE (Ω = √(g/R))"
            ))

            stability_eq0: str
            if all(v is not None for v in [R_val, Omega_val, g_val]):
                coeff_num = float(coeff_eq0_simplified)
                if coeff_num > 0:
                    stability_eq0 = "STABLE"
                    omega_n_eq0 = float(sp.sqrt(coeff_eq0_simplified))
                    final_answer["omega_n_theta0"] = round(omega_n_eq0, 8)
                elif coeff_num < 0:
                    stability_eq0 = "UNSTABLE"
                else:
                    stability_eq0 = "MARGINALLY STABLE"
                steps.append((
                    "Stability at θ=0 (numeric)",
                    f"κ₀ = {coeff_num:.6f}  →  {stability_eq0}"
                ))
            else:
                stability_eq0 = "STABLE if Ω < √(g/R); UNSTABLE if Ω > √(g/R)"
                steps.append((
                    "Critical speed for θ=0 stability",
                    f"Ω_crit = √(g/R)  →  STABLE for Ω < Ω_crit"
                ))
            final_answer["stability_theta0"] = stability_eq0
            final_answer["linearised_coefficient_theta0"] = str(coeff_eq0_simplified)

        # Step 5 – Linearise about non-trivial equilibrium
        if eq_choice in ("all", "theta=arccos"):
            steps.append((
                "Linearisation about θ_eq = arccos(g/(R·Ω²))",
                "Let θ = θ_eq + ε\n"
                "Expand EOM to first order in ε:\n"
                "  θ'' → ε''\n"
                "  sin(θ_eq+ε)·[(g/R) − Ω²cos(θ_eq+ε)]\n"
                "  ≈ [sin(θ_eq) + ε·cos(θ_eq)]·[(g/R) − Ω²cos(θ_eq) + ε·Ω²sin(θ_eq)]\n"
                "At equilibrium: (g/R) − Ω²·cos(θ_eq) = 0\n"
                "Collecting ε terms:\n"
                "  ε'' + Ω²·sin²(θ_eq)·ε = 0"
            ))

            coeff_eq2 = Omega_sym**2 * (1 - cos_eq2**2)   # Omega^2 * sin^2(theta_eq)
            coeff_eq2_simplified = sp.simplify(coeff_eq2)
            omega_n_eq2_sym = sp.sqrt(Omega_sym**2 * (1 - cos_eq2**2))

            steps.append((
                "Natural frequency at non-trivial equilibrium",
                f"ω_n = Ω·|sin(θ_eq)|\n"
                f"    = Ω·√(1 − [g/(R·Ω²)]²)\n"
                f"    = {sp.sstr(sp.simplify(omega_n_eq2_sym))}\n"
                f"Restoring coefficient: Ω²·sin²(θ_eq) = {sp.sstr(coeff_eq2_simplified)}"
            ))

            if all(v is not None for v in [R_val, Omega_val, g_val]):
                cos_val = float(cos_eq2_simplified)
                if abs(cos_val) <= 1.0:
                    import math
                    sin2_val = 1.0 - cos_val**2
                    coeff_num2 = float(Omega_val)**2 * sin2_val  # type: ignore[arg-type]
                    omega_n_num2 = math.sqrt(coeff_num2)
                    stability_eq2 = "STABLE" if coeff_num2 > 0 else "UNSTABLE"
                    steps.append((
                        "Stability at non-trivial equilibrium (numeric)",
                        f"sin²(θ_eq) = 1 − {cos_val:.6f}² = {sin2_val:.6f}\n"
                        f"κ = Ω²·sin²(θ_eq) = {float(Omega_val)**2:.4f} × {sin2_val:.6f} = {coeff_num2:.6f}\n"  # type: ignore[arg-type]
                        f"ω_n = {omega_n_num2:.6f} rad/s\n"
                        f"→ {stability_eq2}"
                    ))
                    final_answer["stability_nontrivial"] = stability_eq2
                    final_answer["omega_n_nontrivial"] = round(omega_n_num2, 8)
                else:
                    steps.append((
                        "Non-trivial equilibrium does not exist",
                        "Ω² < g/R; skip this equilibrium."
                    ))
                    final_answer["stability_nontrivial"] = "N/A (equilibrium does not exist)"
            else:
                final_answer["stability_nontrivial"] = "STABLE (for Ω > √(g/R))"
                final_answer["omega_n_nontrivial"] = "Omega*sqrt(1-(g/(R*Omega^2))^2)"

            final_answer["linearised_coefficient_nontrivial"] = str(coeff_eq2_simplified)

        # Step 6 – Summary sanity check
        sanity_parts: list[str] = []
        sanity_parts.append(
            "EOM verified: θ=0 satisfies sin(0)=0 → EOM vanishes ✓"
        )
        if all(v is not None for v in [R_val, Omega_val, g_val]):
            import math
            omega_crit = math.sqrt(float(g_val) / float(R_val))  # type: ignore[arg-type]
            sanity_parts.append(
                f"Critical speed Ω_crit = √(g/R) = √({float(g_val)}/{float(R_val)}) "  # type: ignore[arg-type]
                f"= {omega_crit:.6f} rad/s"
            )
            if float(Omega_val) > omega_crit:  # type: ignore[arg-type]
                sanity_parts.append(
                    f"Ω = {float(Omega_val):.4f} > Ω_crit → θ=0 UNSTABLE, non-trivial eq. STABLE ✓"  # type: ignore[arg-type]
                )
            else:
                sanity_parts.append(
                    f"Ω = {float(Omega_val):.4f} < Ω_crit → θ=0 STABLE, non-trivial eq. does not exist ✓"  # type: ignore[arg-type]
                )

        return SolverResult(
            problem_type="Rotating Hoop – Lagrange / Equilibrium / Stability",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    def get_input_template(self) -> dict:
        """Return the expected input parameter schema."""
        return {
            "R": {
                "type": "float | None",
                "description": "Hoop radius (m). None → symbolic R.",
                "example": 0.5,
            },
            "m": {
                "type": "float | None",
                "description": "Bead mass (kg). None → symbolic m.",
                "example": 1.0,
            },
            "Omega": {
                "type": "float | None",
                "description": "Hoop angular velocity (rad/s). None → symbolic Omega.",
                "example": 7.0,
            },
            "g": {
                "type": "float | None",
                "description": "Gravitational acceleration (m/s²). None → symbolic g.",
                "example": 9.81,
            },
            "equilibrium_to_analyze": {
                "type": "str",
                "description": (
                    '"all"          → analyse both equilibria (default)\n'
                    '"theta=0"      → analyse only the trivial equilibrium\n'
                    '"theta=arccos" → analyse only the non-trivial equilibrium'
                ),
                "default": "all",
            },
        }
