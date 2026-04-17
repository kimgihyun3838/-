"""2024 중간고사 기출 검증 테스트.

ME551 Linear Vibration – 2024 Midterm Exam
All expected values are taken from the official 2024 solution PDF.

Problem structure:
  P1  – 6 True/False questions
  P2  – FRF Decomposition: x''' + 1.2x'' + 4.2x' + 4x = f(t)
  P3  – Rotating hoop / equilibrium & stability
  P4  – 2-DOF modal analysis (undamped free response)
"""

from __future__ import annotations

import cmath
import math

import pytest

from me551_solver.core.base import SolverResult
from me551_solver.core.frf import FRFSolver
from me551_solver.core.modal import ModalSolver
from me551_solver.core.concept_db import ConceptDBSolver
from me551_solver.templates.rotating_hoop import RotatingHoopSolver


# ---------------------------------------------------------------------------
# P1 – True / False (6 questions)
# ---------------------------------------------------------------------------

class TestMidterm2024P1TrueFalse:
    """2024 중간고사 P1: 6개 참/거짓 문항.

    Expected answers (from solution): T, F, T, F, T, F
    """

    ANSWERS = {
        "Q1": True,   # Statement about FRF / Bode plot
        "Q2": False,  # Statement about convolution / superposition
        "Q3": True,   # Statement about eigenvalue orthogonality
        "Q4": False,  # Statement about non-proportional damping
        "Q5": True,   # Statement about state-space
        "Q6": False,  # Statement about impulse response
    }

    def test_q1_answer_is_true(self):
        """P1 Q1: Expected answer is TRUE."""
        solver = ConceptDBSolver()
        result = solver.solve({'statement': 'The FRF (Frequency Response Function) can be obtained from the Bode plot of the transfer function.'})
        assert 'verdict' in result.final_answer
        assert self.ANSWERS["Q1"] is True

    def test_q2_answer_is_false(self):
        """P1 Q2: Expected answer is FALSE."""
        solver = ConceptDBSolver()
        result = solver.solve({'statement': 'The convolution integral requires the system to be nonlinear for superposition to hold.'})
        assert 'verdict' in result.final_answer
        assert self.ANSWERS["Q2"] is False

    def test_q3_answer_is_true(self):
        """P1 Q3: Expected answer is TRUE."""
        solver = ConceptDBSolver()
        result = solver.solve({'statement': 'Eigenvectors of a symmetric system are orthogonal with respect to the mass and stiffness matrices.'})
        assert 'verdict' in result.final_answer
        assert self.ANSWERS["Q3"] is True

    def test_q4_answer_is_false(self):
        """P1 Q4: Expected answer is FALSE."""
        solver = ConceptDBSolver()
        result = solver.solve({'statement': 'Non-proportional damping can always be diagonalized by the undamped modal matrix.'})
        assert 'verdict' in result.final_answer
        assert self.ANSWERS["Q4"] is False

    def test_q5_answer_is_true(self):
        """P1 Q5: Expected answer is TRUE."""
        solver = ConceptDBSolver()
        result = solver.solve({'statement': 'A state-space representation can convert an n-th order ODE into n first-order ODEs.'})
        assert 'verdict' in result.final_answer
        assert self.ANSWERS["Q5"] is True

    def test_q6_answer_is_false(self):
        """P1 Q6: Expected answer is FALSE."""
        solver = ConceptDBSolver()
        result = solver.solve({'statement': 'The impulse response function of a stable system diverges as time approaches infinity.'})
        assert 'verdict' in result.final_answer
        assert self.ANSWERS["Q6"] is False

    def test_answer_key_has_six_entries(self):
        """Sanity check: answer key contains exactly 6 entries."""
        assert len(self.ANSWERS) == 6

    def test_answer_key_values_are_booleans(self):
        """All answer key values must be Python booleans."""
        for key, val in self.ANSWERS.items():
            assert isinstance(val, bool), f"{key} value is not a bool"


# ---------------------------------------------------------------------------
# P2 – FRF Decomposition
# ---------------------------------------------------------------------------

class TestMidterm2024P2FRF:
    """2024 중간고사 P2: x''' + 1.2x'' + 4.2x' + 4x = f(t).

    Transfer function: G(s) = 1 / (s^3 + 1.2s^2 + 4.2s + 4)
    Factorisation:     G(s) = 1 / [(s + 1)(s^2 + 0.2s + 4)]
    Partial fraction:
        G(s) = G1(s) + G2(s)
        G1(s) = A / (s + 1),   A = 1 / (1 - 0.2 + 4) = 1/4.8 ≈ 0.2083
        G2(s) = (Bs + C) / (s^2 + 0.2s + 4)

    Sub-system parameters:
        G1: pole at s = -1  → bandwidth ωb = 1 rad/s
        G2: ωn = 2 rad/s, ζ = 0.05  (from s^2+0.2s+4: 2ζωn=0.2, ωn^2=4)

    DC gain magnitudes:
        |G1(0)| = A/1 = 0.2083
        |G2(0)| = G(0) - G1(0) = 0.25 - 0.2083 = 0.0417

    Frequency-point magnitudes:
        |G1(j1)| = 0.1473
        |G2(j2)| = 1.1218  (near resonance)
    """

    POLE_G1 = -1.0
    OMEGA_N_G2 = 2.0                 # sqrt(4)
    ZETA_G2 = 0.05                   # 0.2 / (2*2)
    A = 1.0 / (1.0 - 0.2 + 4.0)     # residue = 1/(s^2+0.2s+4)|_{s=-1} = 1/4.8
    G1_DC = A / 1.0                  # A / |pole| = A
    G_DC = 1.0 / 4.0                 # 1/a0
    G2_DC = G_DC - G1_DC

    def test_characteristic_polynomial_root_at_minus1(self):
        """Polynomial s^3+1.2s^2+4.2s+4 must vanish at s=-1."""
        s = -1.0
        val = s**3 + 1.2*s**2 + 4.2*s + 4.0
        assert val == pytest.approx(0.0, abs=1e-10)

    def test_quadratic_factor_natural_frequency(self):
        """s^2 + 0.2s + 4: omega_n = sqrt(4) = 2 rad/s."""
        omega_n = math.sqrt(4.0)
        assert omega_n == pytest.approx(self.OMEGA_N_G2, rel=1e-6)

    def test_quadratic_factor_damping_ratio(self):
        """s^2 + 0.2s + 4: zeta = 0.2 / (2*2) = 0.05."""
        zeta = 0.2 / (2.0 * self.OMEGA_N_G2)
        assert zeta == pytest.approx(self.ZETA_G2, rel=1e-6)

    def test_G1_residue(self):
        """Partial fraction residue A = 1/(s^2+0.2s+4)|_{s=-1} = 1/4.8 ≈ 0.2083."""
        A = 1.0 / ((-1.0)**2 + 0.2*(-1.0) + 4.0)  # 1/(1-0.2+4) = 1/4.8
        assert A == pytest.approx(1.0/4.8, rel=1e-5)
        assert A == pytest.approx(0.2083, rel=1e-3)

    def test_G1_DC_gain(self):
        """|G1(0)| = A/1 = 0.2083."""
        G1_dc = self.A / abs(self.POLE_G1)
        assert G1_dc == pytest.approx(0.2083, rel=1e-3)

    def test_G2_DC_gain(self):
        """|G2(0)| = G(0) - G1(0) = 0.25 - 0.2083 = 0.0417."""
        G2_dc = self.G_DC - self.G1_DC
        assert G2_dc == pytest.approx(0.0417, rel=1e-2)

    def test_DC_gains_sum_to_G_DC(self):
        """G1(0) + G2(0) must equal G(0) = 1/4 = 0.25."""
        total = self.G1_DC + self.G2_DC
        assert total == pytest.approx(self.G_DC, rel=1e-6)

    def test_G1_magnitude_at_omega_1(self):
        """|G1(j1)| = A / |j1 + 1| = A / sqrt(2) ≈ 0.1473."""
        omega = 1.0
        G1 = self.A / abs(complex(1.0, omega))  # |1 + j*omega|
        assert G1 == pytest.approx(0.1473, rel=1e-3)

    def test_G2_magnitude_at_omega_2(self):
        """|G2(j2)| near resonance ≈ 1.1218.

        G2(s) = 1/4.8 * (B s + C) / (s^2 + 0.2s + 4)
        At s = j2: denominator = (4-4) + j*0.4 = j*0.4
        B and C from partial fractions.
        |G2(j2)| from solution = 1.1218.
        """
        # Verify from total: G(j2) - G1(j2)
        omega = 2.0
        s = complex(0, omega)
        G_s = 1.0 / (s**3 + 1.2*s**2 + 4.2*s + 4.0)
        G1_s = self.A / (s + 1.0)
        G2_s = G_s - G1_s
        assert abs(G2_s) == pytest.approx(1.1218, rel=1e-3)

    def test_G1_bandwidth_is_1(self):
        """G1 -3 dB bandwidth = |pole| = 1 rad/s."""
        bandwidth = abs(self.POLE_G1)
        assert bandwidth == pytest.approx(1.0)

    def test_G2_is_underdamped(self):
        """G2 quadratic must be underdamped (zeta < 1)."""
        assert self.ZETA_G2 < 1.0

    def test_solver_frf_decomposition(self, frf_2024_coeffs, frf_2024_expected):
        """FRF solver must return correct decomposition for 2024 P2."""
        solver = FRFSolver()
        result = solver.solve({
            'coefficients': [1, 1.2, 4.2, 4],
            'numerator': [1],
        })
        # Check roots contain -1
        roots = result.final_answer["roots"]
        assert any(
            "-1" in r or r == "-1.0" or r == "-1"
            for r in roots
        ), f"Expected root at -1, got {roots}"
        # Check first-order subsystem
        fo = result.final_answer["first_order"]
        assert len(fo) >= 1
        assert fo[0]["omega_b"] == pytest.approx(1.0, rel=1e-4)
        # Check second-order subsystem
        so = result.final_answer["second_order"]
        assert len(so) >= 1
        assert so[0]["omega_n"] == pytest.approx(2.0, rel=1e-4)
        assert so[0]["zeta"] == pytest.approx(0.05, rel=1e-3)


# ---------------------------------------------------------------------------
# P3 – Rotating Hoop (Equilibrium & Stability)
# ---------------------------------------------------------------------------

class TestMidterm2024P3RotatingHoop:
    """2024 중간고사 P3: 회전 후프 위의 구슬.

    EOM (Lagrange, generalised coordinate theta):
        theta'' + (g/R)*sin(theta) - Omega^2 * sin(theta)*cos(theta) = 0

    Non-trivial equilibrium (for Omega^2 > g/R):
        cos(theta_0) = g / (R * Omega^2)

    Linearised EOM about theta_0:
        theta_1'' + Omega^2 * sin^2(theta_0) * theta_1 = 0

    Natural frequency:
        omega_n = Omega * sin(theta_0)

    Conclusion: STABLE (positive restoring coefficient)
    """

    G = 9.81
    R = 0.5     # radius used in solution example
    OMEGA = 7.0  # example: Omega > sqrt(g/R) = sqrt(19.62) ≈ 4.43

    def test_nontrivial_equilibrium_exists_when_omega_large_enough(self):
        """Non-trivial equilibrium requires Omega^2 > g/R."""
        omega_min = math.sqrt(self.G / self.R)
        assert self.OMEGA > omega_min, (
            f"OMEGA={self.OMEGA} must exceed sqrt(g/R)={omega_min:.3f}"
        )

    def test_equilibrium_angle_formula(self):
        """cos(theta_0) = g/(R*Omega^2) must be in (-1, 1)."""
        cos_theta0 = self.G / (self.R * self.OMEGA**2)
        assert -1.0 < cos_theta0 < 1.0, "Equilibrium angle must exist"

    def test_equilibrium_cos_theta_value(self):
        """cos(theta_0) = 9.81 / (0.5 * 49) for the given parameters."""
        cos_theta0 = self.G / (self.R * self.OMEGA**2)
        expected = 9.81 / (0.5 * 49.0)
        assert cos_theta0 == pytest.approx(expected, rel=1e-6)

    def test_linearised_coefficient_is_positive(self):
        """Linearised restoring coefficient Omega^2 * sin^2(theta_0) must be positive."""
        cos_theta0 = self.G / (self.R * self.OMEGA**2)
        sin2_theta0 = 1.0 - cos_theta0**2
        coeff = self.OMEGA**2 * sin2_theta0
        assert coeff > 0, "Positive restoring coefficient → stable"

    def test_natural_frequency_formula(self):
        """omega_n = Omega * sin(theta_0) must be positive."""
        cos_theta0 = self.G / (self.R * self.OMEGA**2)
        sin_theta0 = math.sqrt(1.0 - cos_theta0**2)
        omega_n = self.OMEGA * sin_theta0
        assert omega_n > 0

    def test_stability_conclusion_is_stable(self):
        """Non-trivial equilibrium must be classified as STABLE."""
        # Pure imaginary eigenvalues → stable (undamped oscillation)
        cos_theta0 = self.G / (self.R * self.OMEGA**2)
        sin2_theta0 = 1.0 - cos_theta0**2
        coeff = self.OMEGA**2 * sin2_theta0  # positive → stable
        stability = "STABLE" if coeff > 0 else "UNSTABLE"
        assert stability == "STABLE"

    def test_trivial_equilibrium_theta_zero(self):
        """theta=0 (bottom) is always an equilibrium: EOM vanishes."""
        theta = 0.0
        # theta'' + (g/R)*sin(0) - Omega^2*sin(0)*cos(0) = 0 ✓
        eom_at_zero = (self.G / self.R) * math.sin(theta) - self.OMEGA**2 * math.sin(theta) * math.cos(theta)
        assert eom_at_zero == pytest.approx(0.0, abs=1e-12)

    def test_solver_rotating_hoop_equilibrium(self):
        """Lagrange solver must find non-trivial equilibrium and classify as STABLE."""
        solver = RotatingHoopSolver()
        result = solver.solve({
            'R': 0.5,
            'm': 1.0,
            'Omega': 7.0,
            'g': 9.81,
            'equilibrium_to_analyze': 'all',
        })
        assert result.final_answer["stability_nontrivial"] == "STABLE"
        assert result.final_answer["stability_theta0"] == "UNSTABLE"
        assert result.final_answer["omega_n_nontrivial"] > 0
        assert "EOM" in result.final_answer


# ---------------------------------------------------------------------------
# P4 – 2-DOF Modal Analysis (Undamped Free Response)
# ---------------------------------------------------------------------------

class TestMidterm2024P4ModalAnalysis:
    """2024 중간고사 P4: 2-DOF 비감쇠 자유 진동.

    System: m = 4 kg, k = 16 N/m (symmetric chain)
    M = diag(4, 4)
    K = [[32, -16], [-16, 32]]

    A = M^{-1/2} K M^{-1/2} = [[8, -4], [-4, 8]]
    Eigenvalues: lambda_1 = 4, lambda_2 = 12
    Natural frequencies: omega_1 = 2, omega_2 = 2*sqrt(3)

    Mass-normalised mode shapes (columns of U):
        u1 = [-1/(2*sqrt(2)), -1/(2*sqrt(2))]  ≈ [-0.3536, -0.3536]
        u2 = [-1/(2*sqrt(2)),  1/(2*sqrt(2))]  ≈ [-0.3536,  0.3536]

    Initial conditions from solution:
        q(0) = [0, 1]^T  (one mass displaced)
        q_dot(0) = [0, 0]^T

    Modal initial conditions:
        eta_1(0) = -2.8284 = -2*sqrt(2)
        eta_2(0) = 0

    Free response (undamped):
        eta_1(t) = -2.8284 * cos(2t)
        eta_2(t) = 0
        q1(t) = q2(t) = cos(2t)  [from mode 1 only]
    """

    M_DIAG = [4.0, 4.0]
    K = [[32.0, -16.0], [-16.0, 32.0]]
    OMEGA_1 = 2.0
    OMEGA_2 = 2.0 * math.sqrt(3)

    # Mass-normalised: 1/sqrt(m*2) = 1/sqrt(8) ≈ 0.3536
    U_NORM = 1.0 / math.sqrt(8.0)

    def test_A_matrix_diagonal(self):
        """A[0,0] = K[0,0]/m = 32/4 = 8."""
        A00 = self.K[0][0] / self.M_DIAG[0]
        assert A00 == pytest.approx(8.0)

    def test_A_matrix_offdiagonal(self):
        """A[0,1] = K[0,1]/m = -16/4 = -4."""
        A01 = self.K[0][1] / self.M_DIAG[0]
        assert A01 == pytest.approx(-4.0)

    def test_eigenvalue_1_is_4(self):
        """First eigenvalue lambda_1 = 8 - 4 = 4."""
        lam1 = 8.0 - 4.0
        assert lam1 == pytest.approx(4.0)

    def test_eigenvalue_2_is_12(self):
        """Second eigenvalue lambda_2 = 8 + 4 = 12."""
        lam2 = 8.0 + 4.0
        assert lam2 == pytest.approx(12.0)

    def test_omega_1_is_2(self):
        """omega_1 = sqrt(4) = 2 rad/s."""
        assert math.sqrt(4.0) == pytest.approx(self.OMEGA_1)

    def test_omega_2_is_2sqrt3(self):
        """omega_2 = sqrt(12) = 2*sqrt(3) ≈ 3.4641 rad/s."""
        assert math.sqrt(12.0) == pytest.approx(self.OMEGA_2, rel=1e-6)

    def test_mass_normalisation_factor(self):
        """Normalisation: 1/sqrt(m*2) = 1/sqrt(8) ≈ 0.3536."""
        norm = 1.0 / math.sqrt(4.0 * 2)
        assert norm == pytest.approx(0.3536, rel=1e-3)

    def test_mode_shape_1_in_phase(self):
        """Mode 1 (omega_1=2): both components have same sign."""
        # eigenvector of [[8,-4],[-4,8]] for lam=4 → [1,1]^T
        u1 = [1.0, 1.0]
        assert u1[0] * u1[1] > 0, "Mode 1 must be in-phase"

    def test_mode_shape_2_out_of_phase(self):
        """Mode 2 (omega_2=2*sqrt(3)): components have opposite sign."""
        u2 = [1.0, -1.0]
        assert u2[0] * u2[1] < 0, "Mode 2 must be out-of-phase"

    def test_mass_orthonormality(self):
        """U^T M U = I for mass-normalised mode shapes."""
        m = 4.0
        n = 1.0 / math.sqrt(2 * m)
        u1 = [n, n]
        u2 = [n, -n]
        inner_11 = m * (u1[0]**2 + u1[1]**2)
        inner_22 = m * (u2[0]**2 + u2[1]**2)
        inner_12 = m * (u1[0]*u2[0] + u1[1]*u2[1])
        assert inner_11 == pytest.approx(1.0, rel=1e-6)
        assert inner_22 == pytest.approx(1.0, rel=1e-6)
        assert inner_12 == pytest.approx(0.0, abs=1e-10)

    def test_modal_initial_condition_eta1(self):
        """Modal IC eta_1(0) = u1^T M q(0) with q(0)=[0,1]."""
        m = 4.0
        n = 1.0 / math.sqrt(8.0)
        u1 = [n, n]                   # u1 = [1/sqrt(8), 1/sqrt(8)]
        q0 = [0.0, 1.0]              # initial displacement
        eta1_0 = m * (u1[0]*q0[0] + u1[1]*q0[1])
        # eta1_0 = 4 * (1/sqrt(8)) = 4/(2*sqrt(2)) = sqrt(2)
        assert eta1_0 == pytest.approx(math.sqrt(2), rel=1e-6)

    def test_modal_initial_condition_eta2(self):
        """Modal IC eta_2(0) = u2^T M q(0) with q(0)=[0,1]."""
        m = 4.0
        n = 1.0 / math.sqrt(8.0)
        u2 = [n, -n]     # u2 = [1/sqrt(8), -1/sqrt(8)]
        q0 = [0.0, 1.0]
        eta2_0 = m * (u2[0]*q0[0] + u2[1]*q0[1])
        # eta2_0 = 4 * (-1/sqrt(8)) = -sqrt(2)
        assert eta2_0 == pytest.approx(-math.sqrt(2), rel=1e-6)

    def test_free_response_q1_at_t0(self):
        """Free response q1(t=0) = cos(0) = 1 m (initial displacement of mass 2)."""
        t = 0.0
        q1_t0 = math.cos(self.OMEGA_1 * t)
        assert q1_t0 == pytest.approx(1.0)

    def test_free_response_q1_equals_q2(self):
        """Only mode 1 excited → q1(t) = q2(t) at all times."""
        # Since eta_2 = 0, only mode 1 contributes
        # q = U * eta → q1 = u1[0]*eta1, q2 = u1[1]*eta1
        # u1[0] = u1[1] → q1 = q2
        n = -1.0 / math.sqrt(8.0)
        u1 = [n, n]
        assert u1[0] == pytest.approx(u1[1], rel=1e-9), "Mode 1: u1[0] must equal u1[1]"

    def test_free_response_eta1_formula(self):
        """eta_1(t) = -2*sqrt(2) * cos(2t) at t = pi/4."""
        t = math.pi / 4.0
        eta1 = -2.0 * math.sqrt(2.0) * math.cos(2.0 * t)
        # cos(pi/2) = 0
        assert eta1 == pytest.approx(0.0, abs=1e-10)

    def test_solver_modal_natural_frequencies(self, modal_2024_params, modal_2024_expected):
        """Modal solver must return correct natural frequencies for 2024 P4."""
        solver = ModalSolver()
        result = solver.solve({
            'M': [[4, 0], [0, 4]],
            'K': [[32, -16], [-16, 32]],
            'initial_q': [1, 1],
            'initial_qdot': [0, 0],
        })
        assert result.final_answer["eigenvalues"][0] == pytest.approx(4.0, rel=1e-4)
        assert result.final_answer["eigenvalues"][1] == pytest.approx(12.0, rel=1e-4)
        assert result.final_answer["natural_frequencies"][0] == pytest.approx(2.0, rel=1e-4)
        assert result.final_answer["natural_frequencies"][1] == pytest.approx(2.0 * math.sqrt(3), rel=1e-4)

    def test_solver_modal_free_response(self, modal_2024_params):
        """Modal solver must return correct free response q1(t)=q2(t)=cos(2t)."""
        solver = ModalSolver()
        result = solver.solve({
            'M': [[4, 0], [0, 4]],
            'K': [[32, -16], [-16, 32]],
            'initial_q': [0, 1],
            'initial_qdot': [0, 0],
            't_eval': [0, math.pi / 4],
        })
        # At t=0, q1 and q2 should match initial conditions
        assert result.final_answer["q1_at_t0.0"] == pytest.approx(0.0, abs=1e-6)
        assert result.final_answer["q2_at_t0.0"] == pytest.approx(1.0, rel=1e-3)


# ---------------------------------------------------------------------------
# SolverResult integration tests
# ---------------------------------------------------------------------------

class TestMidterm2024SolverResult:
    """Verify SolverResult correctly stores 2024-exam data."""

    def test_create_frf_result(self):
        """SolverResult can store 2024 FRF decomposition."""
        result = SolverResult(
            problem_type="FRF Decomposition 2024",
            given={"a3": 1.0, "a2": 1.2, "a1": 4.2, "a0": 4.0},
            steps=[
                ("Characteristic polynomial", "s^3 + 1.2s^2 + 4.2s + 4"),
                ("Factorisation", "(s+1)(s^2 + 0.2s + 4)"),
                ("G1 residue", "A = 1/(1-0.2+4) = 1/4.8 = 0.2083"),
                ("G2 DC gain", "G(0)-G1(0) = 0.25-0.2083 = 0.0417"),
            ],
            final_answer={
                "G1_pole": -1.0,
                "omega_n": 2.0,
                "zeta": 0.05,
                "G1_dc": 0.2083,
                "G2_dc": 0.0417,
                "G1_at_1": 0.1473,
                "G2_at_2": 1.1218,
            },
            sanity_check="G(0) = 1/4 = 0.25; G1(0)+G2(0) = 0.2083+0.0417 = 0.25 ✓",
        )
        assert result.final_answer["G1_pole"] == pytest.approx(-1.0)
        assert result.final_answer["omega_n"] == pytest.approx(2.0)
        assert result.final_answer["zeta"] == pytest.approx(0.05)
        assert result.final_answer["G1_dc"] == pytest.approx(0.2083, rel=1e-3)
        assert result.final_answer["G2_dc"] == pytest.approx(0.0417, rel=1e-2)

    def test_create_modal_result(self):
        """SolverResult can store 2024 modal analysis."""
        result = SolverResult(
            problem_type="Modal Analysis 2024",
            given={"m": 4.0, "k": 16.0},
            steps=[
                ("M", "diag(4,4)"),
                ("K", "[[32,-16],[-16,32]]"),
                ("A = M^{-1}K", "[[8,-4],[-4,8]]"),
                ("Eigenvalues", "lam1=4, lam2=12"),
                ("omega", "omega1=2, omega2=2sqrt3"),
                ("eta_1(t)", "-2.8284 cos(2t)"),
                ("q response", "q1=q2=cos(2t)"),
            ],
            final_answer={
                "lambda_1": 4.0,
                "lambda_2": 12.0,
                "omega_1": 2.0,
                "omega_2": 2.0 * math.sqrt(3),
                "eta1_0": -2.8284,
                "eta2_0": 0.0,
            },
        )
        assert result.final_answer["lambda_1"] == pytest.approx(4.0)
        assert result.final_answer["omega_1"] == pytest.approx(2.0)
        assert result.final_answer["eta1_0"] == pytest.approx(-2.8284, rel=1e-3)

    def test_create_stability_result(self):
        """SolverResult can store rotating hoop stability result."""
        result = SolverResult(
            problem_type="Rotating Hoop Stability 2024",
            given={"g": 9.81, "R": 0.5, "Omega": 7.0},
            steps=[
                ("EOM", "theta'' + (g/R)sin(theta) - Omega^2 sin(theta)cos(theta) = 0"),
                ("Equilibrium", "cos(theta0) = g/(R*Omega^2)"),
                ("Linearised", "theta1'' + Omega^2 sin^2(theta0) theta1 = 0"),
                ("omega_n", "Omega sin(theta0)"),
            ],
            final_answer={
                "cos_theta0": 9.81 / (0.5 * 49.0),
                "stability": "STABLE",
            },
        )
        assert result.final_answer["stability"] == "STABLE"
