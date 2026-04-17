"""Router classification tests and end-to-end solver pipeline tests."""

import pytest
from me551_solver.core.router import classify, get_menu_number


# ======================================================================
# Router classification tests
# ======================================================================

class TestRouter:
    """Test problem type classification."""

    @pytest.mark.parametrize("text,expected", [
        ("x''' + 3x'' + 6x' + 8x = f(t)", "frf"),
        ("Transfer function G(s) partial fraction decomposition", "frf"),
        ("bode plot frequency response ODE", "frf"),
        ("T = 1/2*m*rdot^2, V = -mgR*cos(theta), equilibrium", "lagrange"),
        ("Lagrange equation linearization stability", "lagrange"),
        ("M = [[9,0],[0,9]], K eigenvalue", "modal"),
        ("modal analysis natural frequency mode shape", "modal"),
        ("C = alpha*M + beta*K damping ratio", "damping"),
        ("Rayleigh damping proportional", "damping"),
        ("gyroscopic flutter divergence", "extended"),
        ("skew-symmetric nonconservative circulatory", "extended"),
        ("Is this true or false", "tf"),
        ("T/F superposition linear system", "tf"),
    ])
    def test_classification(self, text, expected):
        result = classify(text)
        assert result.route == expected, \
            f"Expected {expected}, got {result.route} for: {text}"

    def test_menu_numbers(self):
        assert get_menu_number("tf") == "1"
        assert get_menu_number("frf") == "2"
        assert get_menu_number("state_space") == "3"
        assert get_menu_number("inverse_laplace") == "4"
        assert get_menu_number("lagrange") == "5"
        assert get_menu_number("modal") == "6"
        assert get_menu_number("damping") == "7"
        assert get_menu_number("extended") == "8"

    def test_unknown_input(self):
        result = classify("")
        assert result.route == "unknown"
        assert result.confidence == 0.0


# ======================================================================
# End-to-end solver tests
# ======================================================================

class TestE2EFRF:
    """End-to-end: FRF solver pipeline."""

    def test_3rd_order_ode(self):
        from me551_solver.core.frf import FRFSolver
        solver = FRFSolver()
        result = solver.solve({
            "coefficients": [1, 3, 6, 8],
            "numerator": [1],
        })
        fa = result.final_answer
        # Check 1st order + 2nd order subsystems exist
        assert "first_order" in fa
        assert "second_order" in fa
        assert len(fa["first_order"]) >= 1
        assert len(fa["second_order"]) >= 1
        # Check DC gain
        assert abs(fa.get("dc_gain", 0) - 0.125) < 0.01


class TestE2EModal:
    """End-to-end: Modal analysis pipeline."""

    def test_2dof_symmetric(self):
        from me551_solver.core.modal import ModalSolver
        solver = ModalSolver()
        result = solver.solve({
            "M": [[9, 0], [0, 9]],
            "K": [[72, -36], [-36, 36]],
        })
        fa = result.final_answer
        assert "natural_frequencies" in fa
        freqs = fa["natural_frequencies"]
        assert len(freqs) == 2
        # Check frequencies are positive and distinct
        assert freqs[0] > 0
        assert freqs[1] > freqs[0]


class TestE2ELagrange:
    """End-to-end: Lagrange N-DOF pipeline."""

    def test_2dof_spring_mass(self):
        from me551_solver.core.lagrange import LagrangeSolver
        solver = LagrangeSolver()
        result = solver.solve({
            "T_expr": "Rational(1,2)*m*q1dot**2 + Rational(1,2)*m*q2dot**2",
            "V_expr": "Rational(1,2)*k*q1**2 + Rational(1,2)*k*(q2-q1)**2",
            "coords": [("q1", "q1dot"), ("q2", "q2dot")],
            "parameters": ["m", "k"],
            "parameter_values": {"m": 1, "k": 1},
            "find": ["eom", "equilibrium", "linearization", "stability"],
        })
        fa = result.final_answer
        assert "eom" in fa
        assert len(fa["eom"]) == 2
        assert "equilibrium_points" in fa
        assert "linearized" in fa
        assert "stability" in fa


class TestE2EExtended:
    """End-to-end: Extended (gyroscopic) pipeline."""

    def test_gyroscopic_stable(self):
        from me551_solver.core.extended import ExtendedSolver
        solver = ExtendedSolver()
        result = solver.solve({
            "M": [[1, 0], [0, 1]],
            "K": [[3, -1], [-1, 3]],
            "G": [[0, 2], [-2, 0]],
        })
        fa = result.final_answer
        assert fa["system_type"] == "pure_gyroscopic"
        assert "stable" in fa["stability"].lower()


class TestE2EConceptDB:
    """End-to-end: ConceptDB v2 pipeline."""

    def test_true_statement(self):
        from me551_solver.core.concept_db import ConceptDBSolver
        solver = ConceptDBSolver()
        result = solver.solve({
            "statement": "A linear system satisfies both homogeneity and additivity."
        })
        assert result.final_answer["adjusted_verdict"] == "TRUE"

    def test_false_with_trap(self):
        from me551_solver.core.concept_db import ConceptDBSolver
        solver = ConceptDBSolver()
        result = solver.solve({
            "statement": "Structural damping can be used for transient response."
        })
        fa = result.final_answer
        # Should detect trap or give FALSE
        assert fa["adjusted_verdict"] == "FALSE" or len(fa.get("conflicts", [])) > 0


class TestE2EReport:
    """End-to-end: Report rendering modes."""

    def test_exam_answer_mode(self):
        from me551_solver.core.report import ReportEngine
        from me551_solver.core.concept_db import ConceptDBSolver

        solver = ConceptDBSolver()
        result = solver.solve({
            "statement": "Gyroscopic forces do no work on the system."
        })

        # Both render modes should work
        exam = ReportEngine.render_exam_style(result)
        assert "Problem:" in exam

        answer = ReportEngine.render_exam_answer(result)
        assert "판정:" in answer
