"""Shared pytest fixtures for ME551 solver tests."""

from __future__ import annotations

import math

import pytest

from me551_solver.core.base import SolverResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def approx_dict(actual: dict, expected: dict, rel: float = 1e-3) -> bool:
    """Return True when every key in expected matches actual within rel tolerance."""
    for k, v in expected.items():
        assert k in actual, f"Key '{k}' not found in result"
        assert actual[k] == pytest.approx(v, rel=rel), (
            f"Key '{k}': expected {v}, got {actual[k]}"
        )
    return True


# ---------------------------------------------------------------------------
# 2025 Midterm Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def frf_2025_coeffs():
    """2025 P2: x''' + 3x'' + 6x' + 8x = f(t).

    Factored: (s+2)(s^2 + s + 4)
    """
    return {"a3": 1.0, "a2": 3.0, "a1": 6.0, "a0": 8.0}


@pytest.fixture
def frf_2025_expected():
    """Expected FRF decomposition for 2025 P2."""
    return {
        "pole_real": -2.0,          # G1 pole
        "omega_n": 2.0,             # G2 natural frequency
        "zeta": 0.25,               # G2 damping ratio
        "G1_dc_gain": 1.0 / 6.0,   # |G1(0)| = 1/(a0/pole_real_magnitude)
        "bandwidth_G1": 2.0,        # G1 -3dB frequency = |pole|
    }


@pytest.fixture
def modal_2025_params():
    """2025 P4: m=9 kg, k=36 N/m, symmetric 2-DOF chain."""
    return {
        "m": 9.0,
        "k": 36.0,
        "M": [[9.0, 0.0], [0.0, 9.0]],
        "K": [[72.0, -36.0], [-36.0, 72.0]],
        "c": 5.4,   # proportional damping parameter used in problem
    }


@pytest.fixture
def modal_2025_expected():
    """Expected modal analysis results for 2025 P4."""
    return {
        "omega_1": 2.0,
        "omega_2": 2.0 * math.sqrt(3),   # ≈ 3.4641
        "zeta_1": 0.075,   # given proportional damping
        "zeta_2": 0.075 * math.sqrt(3),  # proportional to omega
        # Mass-normalised mode shapes (columns of U)
        # u1 = [-1, -1]/sqrt(18), u2 = [-1, 1]/sqrt(18)
        "u1_ratio": 1.0,    # |u1[0]/u1[1]| = 1
        "u2_ratio": -1.0,   # u2[0]/u2[1] = -1
    }


# ---------------------------------------------------------------------------
# 2024 Midterm Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def frf_2024_coeffs():
    """2024 P2: x''' + 1.2x'' + 4.2x' + 4x = f(t).

    Factored: (s+1)(s^2 + 0.2s + 4)
    """
    return {"a3": 1.0, "a2": 1.2, "a1": 4.2, "a0": 4.0}


@pytest.fixture
def frf_2024_expected():
    """Expected FRF decomposition for 2024 P2."""
    return {
        "pole_real": -1.0,          # G1 pole
        "omega_n": 2.0,             # G2 natural frequency  sqrt(4)
        "zeta": 0.05,               # G2 damping ratio  0.2/(2*2)
        "G1_dc_gain": pytest.approx(0.2083, rel=1e-3),
        "G2_dc_gain": pytest.approx(0.0417, rel=1e-3),
        "G1_at_omega1": pytest.approx(0.1473, rel=1e-3),  # |G1(j1)|
        "G2_at_omega2": pytest.approx(1.1218, rel=1e-3),  # |G2(j2)|
    }


@pytest.fixture
def modal_2024_params():
    """2024 P4: m=4 kg, k=16 N/m, symmetric 2-DOF chain."""
    return {
        "m": 4.0,
        "k": 16.0,
        "M": [[4.0, 0.0], [0.0, 4.0]],
        "K": [[32.0, -16.0], [-16.0, 32.0]],
    }


@pytest.fixture
def modal_2024_expected():
    """Expected modal analysis results for 2024 P4."""
    return {
        "lambda_1": 4.0,    # eigenvalue
        "lambda_2": 12.0,   # eigenvalue
        "omega_1": 2.0,     # sqrt(4)
        "omega_2": 2.0 * math.sqrt(3),  # sqrt(12)
        # Mass-normalised mode shapes
        "u1_ratio": 1.0,    # |u1[0]/u1[1]| = 1 (in-phase)
        "u2_ratio": -1.0,   # u2[0]/u2[1] = -1 (out-of-phase)
    }


# ---------------------------------------------------------------------------
# Previous Midterms (2020-2023) Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def frf_prev_coeffs():
    """Previous midterm P1: same ODE as 2024 – x''' + 1.2x'' + 4.2x' + 4x = f(t)."""
    return {"a3": 1.0, "a2": 1.2, "a1": 4.2, "a0": 4.0}


@pytest.fixture
def modal_prev_proportional():
    """Previous midterm P4 proportional damping case.

    m1=1, m2=4 kg; k1=25, k2=100 N/m; c1=0.5, c2=2 N·s/m
    """
    return {
        "M": [[1.0, 0.0], [0.0, 4.0]],
        "K": [[125.0, -100.0], [-100.0, 100.0]],
        "C_prop": [[2.5, -2.0], [-2.0, 2.0]],   # proportional c1=0.5, c2=2
        "omega_1": pytest.approx(2.0711, rel=1e-3),
        "omega_2": pytest.approx(12.0711, rel=1e-3),
        "zeta_1": pytest.approx(0.0207, rel=1e-2),
        "zeta_2": pytest.approx(0.1207, rel=1e-2),
    }


@pytest.fixture
def modal_prev_nonproportional():
    """Previous midterm P4 nonproportional damping case.

    m1=1, m2=4 kg; k1=25, k2=100 N/m; c1=2, c2=0.5 N·s/m
    """
    return {
        "omega_1": pytest.approx(2.0743, rel=1e-3),
        "zeta_1": pytest.approx(0.07154, rel=1e-2),
        "omega_2": pytest.approx(12.0522, rel=1e-3),
        "zeta_2": pytest.approx(0.09659, rel=1e-2),
    }


# ---------------------------------------------------------------------------
# Homework Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def hw1_q3_state_space():
    """HW1 Q3: x''' + 2x'' + 4x' + 8x = delta(t).

    State space: A = [[0,1,0],[0,0,1],[-8,-4,-2]], b = [0,0,1]^T
    """
    return {
        "A": [[0, 1, 0], [0, 0, 1], [-8, -4, -2]],
        "b": [0, 0, 1],
        "poles": [-2.0, complex(-0.0, 2.0), complex(-0.0, -2.0)],  # s=-2, s=±2j
    }


@pytest.fixture
def hw1_q4_step_response():
    """HW1 Q4: x'' + 4x' + 9x = u(t).

    Steady state x_ss = 1/9 ≈ 0.1111
    Underdamped: omega_n=3, zeta=2/3, omega_d=sqrt(5)
    """
    return {
        "omega_n": 3.0,
        "zeta": pytest.approx(2.0 / 3.0, rel=1e-4),
        "omega_d": pytest.approx(math.sqrt(5), rel=1e-4),
        "x_ss": pytest.approx(1.0 / 9.0, rel=1e-4),
    }


# ---------------------------------------------------------------------------
# Minimal SolverResult factory (for report tests)
# ---------------------------------------------------------------------------

@pytest.fixture
def minimal_result():
    """A minimal SolverResult for report rendering tests."""
    return SolverResult(
        problem_type="Test Problem",
        given={"m": 1.0, "k": 100.0},
        steps=[("EOM", "m x'' + k x = 0"), ("omega_n", "sqrt(k/m) = 10")],
        final_answer={"omega_n": 10.0},
        sanity_check="omega_n = sqrt(100/1) = 10 rad/s",
    )


@pytest.fixture
def full_result():
    """A fully-populated SolverResult for report tests."""
    return SolverResult(
        problem_type="FRF Decomposition",
        given={"a3": 1.0, "a2": 1.2, "a1": 4.2, "a0": 4.0},
        steps=[
            ("Characteristic polynomial", "s^3 + 1.2s^2 + 4.2s + 4"),
            ("Factorisation", "(s+1)(s^2 + 0.2s + 4)"),
            ("G1", "1/(s+1)  x  (1/4.8)"),
            ("G2", "1/(s^2 + 0.2s + 4)  x  (1/4.8)"),
        ],
        final_answer={
            "G1_pole": -1.0,
            "omega_n": 2.0,
            "zeta": 0.05,
        },
        sanity_check="G(0) = 1/4 = 0.25. G1(0)+G2(0) = 0.2083+0.0417 = 0.25 ✓",
    )
