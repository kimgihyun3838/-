"""Core solver modules."""

from .base import SolverResult, BaseSolver
from .fourier import FourierSolver
from .convolution import ConvolutionSolver

__all__ = ["SolverResult", "BaseSolver", "FourierSolver", "ConvolutionSolver"]
