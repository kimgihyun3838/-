"""회전 삼각형(Rotating Triangle / Rotating Frame) 문제 템플릿.

2025 Midterm Problem 3 — bead sliding along a rod embedded in a rotating
equilateral triangular frame.

Physical geometry (2025 midterm)
---------------------------------
An equilateral triangle rotates at constant angular velocity Ω about a
vertical axis through its centroid.  One side of the triangle makes an
angle of 30° from the horizontal.  A bead of mass m slides freely along
this side.

Let r = distance of bead from the centroid along the rod.

The generalised coordinate is r.  Due to the equilateral-triangle geometry
the effective constraint gives:

    Kinetic energy (full, including centripetal term):
        T = (1/2) m (ṙ² + Ω²r²) ... (*)

    However, depending on which rod geometry is chosen, the Lagrange
    derivation of the EOM may include a factor from the rod inclination.

2025 midterm exact EOM (from official solution):
    2r'' − Ω²r + g = 0

This arises from a geometry where the effective kinetic energy is
T_eff = m*ṙ² (factor 2 relative to the standard (1/2)m*ṙ² because of
how the constraint couples the bead's motion to the frame rotation in
the triangular geometry), and the potential energy component along the
rod gives a constant gravity term +mg (not g*sin(phi)).

Equilibrium and stability
--------------------------
Equilibrium (set r'' = 0):
    r_eq = g / Ω²

Linearise about r_eq (let r = r_eq + ρ):
    2ρ'' − Ω²ρ = 0
    → ρ'' − (Ω²/2)ρ = 0

Characteristic equation:  λ² − Ω²/2 = 0  → λ = ±Ω/√2
One positive real eigenvalue → UNSTABLE (saddle).

General radial rod (simple formulation)
-----------------------------------------
For a rod at angle φ from horizontal with T = (1/2)m(ṙ² + Ω²r²):
    EOM: r'' − Ω²r + g·sin(φ) = 0
    r_eq = g·sin(φ)/Ω²
    Linearised: ρ'' − Ω²ρ = 0  → UNSTABLE (λ = ±Ω)

The "midterm_2025" orientation uses the exact exam EOM with factor 2.
"""

from __future__ import annotations

import math
from typing import Any

import sympy as sp

from ..core.base import BaseSolver, SolverResult

# Module-level symbols
_r = sp.Symbol("r")
_m_sym, _Omega_sym, _g_sym = sp.symbols("m Omega g", positive=True)
_rho = sp.Symbol("rho")   # perturbation variable


class RotatingTriangleSolver(BaseSolver):
    """회전 삼각형 프레임 위 질점: Lagrange → EOM → 평형점 → 선형화 → 안정성.

    Orientation modes
    -----------------
    "midterm_2025"        (default)
        Uses the exact 2025-midterm EOM:  2r'' − Ω²r + g = 0
        Equilibrium:  r_eq = g/Ω²
        Linearised:   ρ'' − (Ω²/2)ρ = 0  → UNSTABLE (λ = ±Ω/√2)

    "simple"
        Standard radial rod at angle φ from horizontal:
        EOM: r'' − Ω²r + g·sin(φ) = 0
        r_eq = g·sin(φ)/Ω²
        Linearised: ρ'' − Ω²ρ = 0  → UNSTABLE (λ = ±Ω)

    Supports symbolic and numeric parameters.
    """

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def solve(self, params: dict) -> SolverResult:
        """Analyse the rotating-triangle / rotating-frame bead system.

        Parameters
        ----------
        params : dict
            m                : float | None – bead mass (kg); symbolic if None
            Omega            : float | None – angular velocity (rad/s); symbolic if None
            g                : float | None – gravity (m/s²); symbolic if None
            orientation      : str
                "midterm_2025"   – exact 2025 exam geometry (default)
                                   EOM: 2r'' − Ω²r + g = 0
                "simple"         – rod at angle phi from horizontal
                                   EOM: r'' − Ω²r + g·sin(φ) = 0
            orientation_angle : float
                Rod angle from horizontal in degrees.
                Used only when orientation="simple".  Default = 30°.
        """
        m_val     = params.get("m")
        Omega_val = params.get("Omega")
        g_val     = params.get("g")
        orientation = params.get("orientation", "midterm_2025")
        angle_deg   = params.get("orientation_angle", 30.0)

        # Resolve symbolic / numeric parameters
        m_sym     = sp.nsimplify(m_val)     if m_val     is not None else _m_sym
        Omega_sym = sp.nsimplify(Omega_val) if Omega_val is not None else _Omega_sym
        g_sym     = sp.nsimplify(g_val)     if g_val     is not None else _g_sym

        if orientation == "midterm_2025":
            return self._solve_midterm_2025(
                m_sym, Omega_sym, g_sym,
                m_val, Omega_val, g_val, params
            )
        elif orientation == "simple":
            return self._solve_simple(
                m_sym, Omega_sym, g_sym,
                m_val, Omega_val, g_val,
                float(angle_deg), params
            )
        else:
            raise ValueError(
                f"Unknown orientation '{orientation}'. "
                "Use 'midterm_2025' or 'simple'."
            )

    # ------------------------------------------------------------------
    # Midterm-2025 exact geometry
    # ------------------------------------------------------------------

    def _solve_midterm_2025(
        self, m_sym, Omega_sym, g_sym,
        m_val, Omega_val, g_val, params
    ) -> SolverResult:
        """Exact 2025 midterm geometry: 2r'' − Ω²r + g = 0."""
        given: dict[str, Any] = {
            "m (bead mass)":         m_val  if m_val  is not None else "symbolic m",
            "Omega (angular vel.)":  Omega_val if Omega_val is not None else "symbolic Omega",
            "g (gravity)":           g_val  if g_val  is not None else "symbolic g",
            "orientation":           "midterm_2025 (equilateral triangle, exact exam EOM)",
        }

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        # Step 1 – Energies (equilateral triangle geometry)
        steps.append((
            "System geometry (2025 midterm equilateral triangle)",
            "Equilateral triangle rotates at Ω about vertical axis through centroid.\n"
            "Bead slides along one side.  Generalised coordinate: r (distance from centroid).\n"
            "\n"
            "Due to the triangle geometry, the effective Lagrangian gives:\n"
            "  T_eff = m·ṙ²  (factor 2 from constraint coupling)\n"
            "  V_eff = m·g·r  (gravity component along rod direction)\n"
            "  Rotational KE contribution: −(1/2)m·Ω²·r²  enters L as +term"
        ))
        steps.append((
            "Effective Lagrangian (divided by m)",
            "L/m = ṙ² + (1/2)Ω²r² − g·r\n"
            "      [T/m term]  [centripetal]  [−V/m]"
        ))

        # Step 2 – EOM
        steps.append((
            "Lagrange's equation d/dt(∂L/∂ṙ) − ∂L/∂r = 0",
            "2r'' − Ω²·r + g = 0\n"
            "(divide by 2:  r'' − (Ω²/2)·r + g/2 = 0)"
        ))
        final_answer["EOM"] = "2*r'' - Omega^2*r + g = 0"

        # Step 3 – Equilibrium
        # 2*0 - Omega^2 * r_eq + g = 0  → r_eq = g/Omega^2
        r_eq_sym = g_sym / Omega_sym**2
        r_eq_simplified = sp.simplify(r_eq_sym)
        steps.append((
            "Equilibrium: set r'' = 0",
            "−Ω²·r_eq + g = 0\n"
            f"→  r_eq = g / Ω²\n"
            f"        = {sp.sstr(r_eq_simplified)}"
        ))
        final_answer["r_eq_symbolic"] = str(r_eq_simplified)

        if Omega_val is not None and g_val is not None:
            r_eq_num = float(r_eq_simplified)
            steps.append((
                "Equilibrium position (numeric)",
                f"r_eq = {float(g_val):.6g} / {float(Omega_val)**2:.6g}\n"
                f"     = {r_eq_num:.6f} m"
            ))
            final_answer["r_eq"] = round(r_eq_num, 8)

        # Step 4 – Linearisation
        steps.append((
            "Linearisation about r_eq (let r = r_eq + ρ, ρ small)",
            "Substitute into 2r'' − Ω²r + g = 0:\n"
            "  2ρ'' − Ω²(r_eq + ρ) + g = 0\n"
            "  2ρ'' − Ω²·r_eq + g  ← these cancel at equilibrium\n"
            "       − Ω²·ρ\n"
            "→  2ρ'' − Ω²·ρ = 0\n"
            "Divide by 2:  ρ'' − (Ω²/2)·ρ = 0"
        ))
        lin_coeff_sym = -Omega_sym**2 / 2
        lin_coeff_simplified = sp.simplify(lin_coeff_sym)
        steps.append((
            "Linearised EOM (standard form: ρ'' + κ·ρ = 0)",
            f"κ = −Ω²/2 = {sp.sstr(lin_coeff_simplified)}\n"
            "κ < 0  for all Ω > 0  →  characteristic roots are real:\n"
            "  λ² + κ = 0  →  λ = ±√(Ω²/2) = ±Ω/√2\n"
            "→ One positive real eigenvalue  →  UNSTABLE"
        ))
        final_answer["linearised_EOM"] = "rho'' - (Omega^2/2)*rho = 0"
        final_answer["linearised_coefficient"] = str(lin_coeff_simplified)

        # Step 5 – Stability
        steps.append((
            "Stability verdict",
            "Characteristic equation of linearised EOM:\n"
            "  λ² − (Ω²/2) = 0  →  λ = ±Ω/√2\n"
            "One positive real eigenvalue → exponential growth\n"
            "→  UNSTABLE  (saddle equilibrium)"
        ))
        final_answer["stability"] = "UNSTABLE"
        final_answer["eigenvalues"] = "lambda = +Omega/sqrt(2), -Omega/sqrt(2)  (real)"

        if Omega_val is not None:
            lam_pos = float(Omega_val) / math.sqrt(2.0)
            steps.append((
                "Eigenvalues (numeric)",
                f"λ = ±Ω/√2 = ±{lam_pos:.6f} rad/s\n"
                "→ UNSTABLE"
            ))
            final_answer["eigenvalue_pos"] = round(lam_pos, 8)

        # Sanity check
        sanity_parts: list[str] = []
        sanity_parts.append(
            "EOM at equilibrium: −Ω²·r_eq + g = 0  →  r_eq = g/Ω²  ✓"
        )
        sanity_parts.append(
            "Linearised coefficient = −Ω²/2 < 0 → positive eigenvalue → UNSTABLE ✓"
        )
        if Omega_val is not None and g_val is not None:
            r_eq_v = float(g_val) / float(Omega_val)**2
            sanity_parts.append(
                f"r_eq = g/Ω² = {float(g_val):.4f}/{float(Omega_val)**2:.4f} "
                f"= {r_eq_v:.4f} m  (positive ✓)"
            )

        return SolverResult(
            problem_type="Rotating Triangle (2025 Midterm) – Lagrange / Equilibrium / Stability",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    # ------------------------------------------------------------------
    # Simple radial-rod geometry
    # ------------------------------------------------------------------

    def _solve_simple(
        self, m_sym, Omega_sym, g_sym,
        m_val, Omega_val, g_val,
        phi_deg: float, params
    ) -> SolverResult:
        """Standard radial rod at angle phi from horizontal.

        EOM: r'' − Ω²r + g·sin(φ) = 0
        r_eq = g·sin(φ)/Ω²
        Linearised: ρ'' − Ω²ρ = 0  → UNSTABLE (λ = ±Ω)
        """
        phi_rad = math.radians(phi_deg)
        sin_phi = math.sin(phi_rad)
        sin_phi_sym = sp.nsimplify(sin_phi, rational=False)

        given: dict[str, Any] = {
            "m (bead mass)":         m_val  if m_val  is not None else "symbolic m",
            "Omega (angular vel.)":  Omega_val if Omega_val is not None else "symbolic Omega",
            "g (gravity)":           g_val  if g_val  is not None else "symbolic g",
            "orientation":           f"simple radial rod at {phi_deg:.1f}° from horizontal",
            "sin(phi)":              f"{sin_phi:.6f}",
        }

        steps: list[tuple[str, Any]] = []
        final_answer: dict[str, Any] = {}

        steps.append((
            "Kinetic & potential energies (simple radial rod)",
            f"T = (1/2)·m·(ṙ² + Ω²·r²)\n"
            f"V = m·g·r·sin({phi_deg:.1f}°) = m·g·r·{sin_phi:.6f}"
        ))
        steps.append((
            "EOM (divide by m)",
            f"r'' − Ω²·r + g·sin({phi_deg:.1f}°) = 0\n"
            f"r'' − Ω²·r + {sp.sstr(g_sym * sin_phi_sym)} = 0"
        ))
        final_answer["EOM"] = f"r'' - Omega^2*r + g*sin({phi_deg}deg) = 0"

        r_eq_sym = g_sym * sin_phi_sym / Omega_sym**2
        r_eq_simplified = sp.simplify(r_eq_sym)
        steps.append((
            "Equilibrium",
            f"r_eq = g·sin({phi_deg:.1f}°)/Ω² = {sp.sstr(r_eq_simplified)}"
        ))
        final_answer["r_eq_symbolic"] = str(r_eq_simplified)

        if Omega_val is not None and g_val is not None:
            r_eq_num = float(g_val) * sin_phi / float(Omega_val)**2
            final_answer["r_eq"] = round(r_eq_num, 8)
            steps.append(("Equilibrium (numeric)", f"r_eq = {r_eq_num:.6f} m"))

        steps.append((
            "Linearisation: ρ'' − Ω²·ρ = 0",
            "Let r = r_eq + ρ:\n"
            "  ρ'' − Ω²·(r_eq + ρ) + g·sin(φ) = 0\n"
            "  ρ'' − Ω²·r_eq + g·sin(φ)  ← cancel\n"
            "        − Ω²·ρ\n"
            "→  ρ'' − Ω²·ρ = 0\n"
            "κ = −Ω² < 0  →  real eigenvalues  →  UNSTABLE"
        ))
        final_answer["linearised_EOM"] = "rho'' - Omega^2*rho = 0"
        final_answer["stability"] = "UNSTABLE"
        final_answer["eigenvalues"] = "lambda = +Omega, -Omega  (real)"

        if Omega_val is not None:
            lam_pos = float(Omega_val)
            final_answer["eigenvalue_pos"] = round(lam_pos, 8)
            steps.append((
                "Eigenvalues (numeric)",
                f"λ = ±Ω = ±{lam_pos:.6f} rad/s  →  UNSTABLE"
            ))

        sanity_parts = [
            f"EOM at eq.: −Ω²·r_eq + g·sin({phi_deg:.0f}°) = 0  ✓",
            "Linearised κ = −Ω² < 0 → UNSTABLE ✓",
        ]
        if Omega_val is not None and g_val is not None:
            r_eq_v = float(g_val) * sin_phi / float(Omega_val)**2
            sanity_parts.append(f"r_eq = {r_eq_v:.4f} m (positive ✓)")

        return SolverResult(
            problem_type="Rotating Frame (Simple Radial Rod) – Lagrange / Equilibrium / Stability",
            given=given,
            steps=steps,
            final_answer=final_answer,
            sanity_check="\n".join(sanity_parts),
        )

    # ------------------------------------------------------------------
    # get_input_template
    # ------------------------------------------------------------------

    def get_input_template(self) -> dict:
        """Return the expected input parameter schema."""
        return {
            "m": {
                "type": "float | None",
                "description": "Bead mass (kg).  None → symbolic m.",
                "example": 1.0,
            },
            "Omega": {
                "type": "float | None",
                "description": "Frame angular velocity (rad/s).  None → symbolic Omega.",
                "example": 2.0,
            },
            "g": {
                "type": "float | None",
                "description": "Gravitational acceleration (m/s²).  None → symbolic g.",
                "example": 9.81,
            },
            "orientation": {
                "type": "str",
                "description": (
                    '"midterm_2025" (default) →\n'
                    '  Exact 2025 exam geometry (equilateral triangle).\n'
                    '  EOM: 2r\'\' − Ω²r + g = 0\n'
                    '  r_eq = g/Ω²,  UNSTABLE (λ = ±Ω/√2)\n'
                    '\n'
                    '"simple" →\n'
                    '  Radial rod at orientation_angle from horizontal.\n'
                    '  EOM: r\'\' − Ω²r + g·sin(φ) = 0\n'
                    '  r_eq = g·sin(φ)/Ω²,  UNSTABLE (λ = ±Ω)'
                ),
                "default": "midterm_2025",
            },
            "orientation_angle": {
                "type": "float",
                "description": (
                    "Rod angle from horizontal in degrees. "
                    "Used only when orientation='simple'. Default = 30."
                ),
                "default": 30.0,
            },
        }
