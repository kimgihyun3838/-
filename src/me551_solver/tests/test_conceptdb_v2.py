"""ConceptDB v2 accuracy tests against 2024+2025 midterm T/F questions.

Tests that the new structural analysis + conflict detection correctly
handles the 12 exam T/F questions (especially the 6 that v1 got wrong).
"""

import pytest
from me551_solver.core.concept_db import ConceptDBSolver


@pytest.fixture(scope="module")
def solver():
    return ConceptDBSolver()


# ======================================================================
# 2024 Midterm T/F (6 questions)
# ======================================================================

class Test2024TF:
    """2024 midterm T/F questions."""

    def test_1a_linear_superposition(self, solver):
        """1a: Linear = superposition, nonlinear = at least one fails → TRUE"""
        result = solver.solve({
            "statement": "A linear system satisfies both homogeneity and additivity. "
                         "A nonlinear system violates at least one of these properties."
        })
        assert result.final_answer["adjusted_verdict"] == "TRUE"

    def test_1b_convolution_both_domains(self, solver):
        """1b: Both time and freq domain use convolution → FALSE
        (freq domain uses multiplication, not convolution)"""
        result = solver.solve({
            "statement": "The response can be obtained by convolution in both time domain "
                         "and frequency domain."
        })
        adj = result.final_answer["adjusted_verdict"]
        conflicts = result.final_answer.get("conflicts", [])
        # Should detect conflict: "both" + "convolution" + "frequency domain"
        assert adj == "FALSE" or len(conflicts) > 0, \
            f"Should detect trap: got verdict={adj}, conflicts={conflicts}"

    def test_1c_virtual_work(self, solver):
        """1c: Virtual work = external force * displacement compatible with constraints → TRUE"""
        result = solver.solve({
            "statement": "In the principle of virtual work, the virtual work is the work done by "
                         "applied forces through virtual displacements compatible with constraints."
        })
        assert result.final_answer["adjusted_verdict"] == "TRUE"

    def test_1d_harmonic_input_nonlinear(self, solver):
        """1d: Harmonic input → harmonic output for both linear and nonlinear → FALSE"""
        result = solver.solve({
            "statement": "For a harmonic input, both linear and nonlinear systems produce "
                         "a harmonic output at the same frequency."
        })
        assert result.final_answer["adjusted_verdict"] == "FALSE"

    def test_1e_lyapunov(self, solver):
        """1e: Negative definite V + positive definite dV/dt → asymptotic stability → TRUE"""
        result = solver.solve({
            "statement": "If a Lyapunov function V is positive definite and its time derivative "
                         "is negative definite, the equilibrium is asymptotically stable."
        })
        assert result.final_answer["adjusted_verdict"] == "TRUE"

    def test_1f_state_dimension_p(self, solver):
        """1f: p-th order n-DOF → p-dimensional state → FALSE (should be p*n)"""
        result = solver.solve({
            "statement": "A p-th order system with n degrees of freedom converted to first-order "
                         "state form has a state vector of dimension p."
        })
        adj = result.final_answer["adjusted_verdict"]
        conflicts = result.final_answer.get("conflicts", [])
        # Should detect numerical mismatch: p vs p*n
        assert adj == "FALSE" or len(conflicts) > 0, \
            f"Should detect dimension trap: got verdict={adj}, conflicts={conflicts}"


# ======================================================================
# 2025 Midterm T/F (6 questions)
# ======================================================================

class Test2025TF:
    """2025 midterm T/F questions."""

    def test_1a_convolution_free_and_forced(self, solver):
        """1a: Convolution/FRF applies to both free and forced vibration → FALSE
        (only forced)"""
        result = solver.solve({
            "statement": "The convolution integral and FRF can be used to compute the response "
                         "for both free and forced vibration of a linear system."
        })
        adj = result.final_answer["adjusted_verdict"]
        conflicts = result.final_answer.get("conflicts", [])
        assert adj == "FALSE" or len(conflicts) > 0, \
            f"Should detect free vibration trap: verdict={adj}, conflicts={conflicts}"

    def test_1b_state_dimension_n(self, solver):
        """1b: 2nd order N-DOF → N-dimensional state → FALSE (should be 2N)"""
        result = solver.solve({
            "statement": "A second-order N-DOF system converted to first-order state-space form "
                         "has an N-dimensional state vector."
        })
        adj = result.final_answer["adjusted_verdict"]
        conflicts = result.final_answer.get("conflicts", [])
        assert adj == "FALSE" or len(conflicts) > 0, \
            f"Should detect dimension trap: verdict={adj}, conflicts={conflicts}"

    def test_1c_virtual_work_internal(self, solver):
        """1c: Virtual work eliminates internal/constraint forces → TRUE"""
        result = solver.solve({
            "statement": "In the principle of virtual work, constraint forces and internal forces "
                         "are eliminated from the equations of motion."
        })
        assert result.final_answer["adjusted_verdict"] == "TRUE"

    def test_1d_structural_damping_transient(self, solver):
        """1d: Structural damping usable for transient → FALSE (harmonic only)"""
        result = solver.solve({
            "statement": "Structural damping can be used for transient response analysis "
                         "as well as harmonic excitation."
        })
        adj = result.final_answer["adjusted_verdict"]
        conflicts = result.final_answer.get("conflicts", [])
        assert adj == "FALSE" or len(conflicts) > 0, \
            f"Should detect transient trap: verdict={adj}, conflicts={conflicts}"

    def test_1e_spd_eigenvalues(self, solver):
        """1e: SPD matrix has positive real eigenvalues, real orthogonal eigenvectors → TRUE"""
        result = solver.solve({
            "statement": "A symmetric positive definite matrix has all positive real eigenvalues "
                         "and its eigenvectors can be chosen to be real and orthogonal."
        })
        assert result.final_answer["adjusted_verdict"] == "TRUE"

    def test_1f_gyroscopic_cannot_diagonalize(self, solver):
        """1f: Gyroscopic system cannot be diagonalized → FALSE
        (can be diagonalized via bi-orthogonality with left/right eigenvectors)"""
        result = solver.solve({
            "statement": "A gyroscopic system cannot be diagonalized using modal analysis."
        })
        adj = result.final_answer["adjusted_verdict"]
        conflicts = result.final_answer.get("conflicts", [])
        assert adj == "FALSE" or len(conflicts) > 0, \
            f"Should detect negation trap: verdict={adj}, conflicts={conflicts}"


# ======================================================================
# Accuracy summary
# ======================================================================

class TestAccuracySummary:
    """Check overall accuracy meets threshold."""

    EXAM_QUESTIONS = [
        # (statement, correct_answer)
        # 2024
        ("A linear system satisfies both homogeneity and additivity.", "TRUE"),
        ("The response can be obtained by convolution in both time domain and frequency domain.", "FALSE"),
        ("In the principle of virtual work, the virtual work is done by applied forces through virtual displacements compatible with constraints.", "TRUE"),
        ("For a harmonic input, both linear and nonlinear systems produce a harmonic output.", "FALSE"),
        ("If a Lyapunov function V is positive definite and dV/dt is negative definite, the equilibrium is asymptotically stable.", "TRUE"),
        ("A p-th order n-DOF system converted to state form has a state vector of dimension p.", "FALSE"),
        # 2025
        ("The convolution integral and FRF can be used for both free and forced vibration.", "FALSE"),
        ("A second-order N-DOF system has an N-dimensional state vector in first-order form.", "FALSE"),
        ("In virtual work, constraint forces and internal forces are eliminated.", "TRUE"),
        ("Structural damping can be used for transient response analysis.", "FALSE"),
        ("A symmetric positive definite matrix has positive real eigenvalues and real orthogonal eigenvectors.", "TRUE"),
        ("A gyroscopic system cannot be diagonalized using modal analysis.", "FALSE"),
    ]

    def test_overall_accuracy(self, solver):
        """Overall accuracy should be >= 75% (9/12)."""
        correct = 0
        total = len(self.EXAM_QUESTIONS)
        details = []

        for stmt, expected in self.EXAM_QUESTIONS:
            result = solver.solve({"statement": stmt})
            adj = result.final_answer.get("adjusted_verdict", "")
            ok = adj == expected
            if ok:
                correct += 1
            details.append(f"{'O' if ok else 'X'} [{expected}] adj={adj} | {stmt[:60]}...")

        accuracy = correct / total
        report = "\n".join(details)
        assert accuracy >= 0.75, (
            f"Accuracy {correct}/{total} = {accuracy:.0%} < 75%\n{report}"
        )

    def test_false_accuracy(self, solver):
        """FALSE question accuracy should be >= 57% (4/7)."""
        false_qs = [(s, e) for s, e in self.EXAM_QUESTIONS if e == "FALSE"]
        correct = 0
        for stmt, expected in false_qs:
            result = solver.solve({"statement": stmt})
            adj = result.final_answer.get("adjusted_verdict", "")
            if adj == expected:
                correct += 1
        accuracy = correct / len(false_qs)
        assert accuracy >= 0.57, (
            f"FALSE accuracy {correct}/{len(false_qs)} = {accuracy:.0%} < 57%"
        )
