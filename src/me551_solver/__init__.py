"""ME551 Linear Vibration Engineering - Midterm Exam Solver

Quick-start
-----------
>>> from me551_solver.templates import RotatingHoopSolver
>>> result = RotatingHoopSolver().solve({"R": 0.5, "Omega": 7.0, "g": 9.81})

>>> from me551_solver.templates import RotatingTriangleSolver
>>> result = RotatingTriangleSolver().solve({"Omega": 2.0, "g": 9.81})

>>> from me551_solver.templates import TwoDOFChainSolver
>>> result = TwoDOFChainSolver().solve({"masses": [9, 9], "springs": [36, 36, 36]})

>>> from me551_solver.templates import FRFDecomposeSolver
>>> result = FRFDecomposeSolver().solve({"ode_string": "x''' + 3x'' + 6x' + 8x = f(t)"})

>>> from me551_solver.core.report import ReportEngine
>>> print(ReportEngine.render_exam_style(result))
"""

__version__ = "0.2.0"
