"""Problem templates – 자주 출제되는 유형의 사전 구성 솔버.

Available templates
-------------------
RotatingHoopSolver
    Bead on a rotating circular hoop (2024 midterm P3).
    Lagrange → EOM → equilibria → linearisation → stability.

RotatingTriangleSolver
    Bead on a rotating triangular/radial frame (2025 midterm P3).
    Lagrange → EOM → equilibrium → linearisation → stability (UNSTABLE).

TwoDOFChainSolver
    Standard 2-DOF spring-mass chain (2024/2025 midterm P4).
    Auto-assembles M and K → eigenvalues → mode shapes → modal response.

FRFDecomposeSolver
    FRF decomposition wrapper (2024/2025 midterm P2).
    Accepts ODE string or coefficient list → delegates to FRFSolver.
"""

from .frf_decompose import FRFDecomposeSolver, parse_ode_string
from .rotating_hoop import RotatingHoopSolver
from .rotating_triangle import RotatingTriangleSolver
from .two_dof_chain import TwoDOFChainSolver

__all__ = [
    "FRFDecomposeSolver",
    "parse_ode_string",
    "RotatingHoopSolver",
    "RotatingTriangleSolver",
    "TwoDOFChainSolver",
]
