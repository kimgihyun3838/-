"""Verify ConceptDBSolver against 60 T/F questions."""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from me551_solver.core.concept_db import ConceptDBSolver

solver = ConceptDBSolver()

questions = [
    (1, "A linear system satisfies homogeneity and additivity, and therefore satisfies the superposition principle.", True),
    (2, "A nonlinear system must violate both homogeneity and additivity simultaneously.", False),
    (3, "For an arbitrary excitation of an LTI system, the response can be obtained by convolution of the input with the impulse response in the time domain.", True),
    (4, "For an arbitrary excitation of an LTI system, the response in the frequency domain is obtained by convolution of the FRF and the input spectrum.", False),
    (5, "The FRF method and the convolution method can be used to obtain the response due to initial conditions only.", False),
    (6, "A second-order linear vibration system with N degrees of freedom can be rewritten as a first-order state equation of dimension 2N.", True),
    (7, "A second-order linear MIMO system with N degrees of freedom can always be rewritten as a first-order state equation of dimension N.", False),
    (8, "For a p-th order linear MIMO system with n degrees of freedom, the first-order state equation generally has dimension pn.", True),
    (9, "The transition matrix of a first-order state equation x_dot=Ax is e^(At).", True),
    (10, "The transition matrix changes the eigenvalues of the system matrix A.", False),
    (11, "The unit step function is the time integral of the Dirac delta function.", True),
    (12, "The unit ramp function is the derivative of the unit step function.", False),
    (13, "For a first-order system, a unit impulse input is mathematically equivalent to an initial displacement of magnitude 1/c.", True),
    (14, "For a second-order system, a unit impulse input is mathematically equivalent to an initial velocity of magnitude 1/m.", True),
    (15, "For a second-order system, a unit impulse input is mathematically equivalent to an initial displacement of magnitude 1/k.", False),
    (16, "A first-order system is typically characterized by a break frequency or time constant.", True),
    (17, "An underdamped second-order system is typically characterized by the undamped natural frequency and the damping ratio.", True),
    (18, "For an underdamped second-order system, the damped natural frequency is larger than the undamped natural frequency.", False),
    (19, "For a harmonic input applied to a linear system, the steady-state output has the same frequency as the input.", True),
    (20, "For a harmonic input applied to a nonlinear system, the steady-state output must have only the input frequency.", False),
    (21, "The compliance FRF is the reciprocal of the dynamic stiffness for an SDOF system.", True),
    (22, "At zero frequency, the compliance of a second-order SDOF system is 1/k.", True),
    (23, "At very high frequency, the compliance magnitude of a second-order SDOF system tends to zero.", True),
    (24, "Increasing damping ratio generally lowers the resonance peak of the compliance FRF.", True),
    (25, "At resonance of an undamped second-order system under harmonic force, the steady-state amplitude remains finite.", False),
    (26, "Structural damping is commonly associated with hysteretic energy loss of elastic materials under cyclic stress.", True),
    (27, "Structural damping is usually represented by an equivalent viscous damping only for harmonic analysis.", True),
    (28, "Structural damping can be used in transient analysis with the same validity as in harmonic analysis.", False),
    (29, "For structural damping, the energy dissipation per cycle is proportional to the square of the response amplitude.", True),
    (30, "The loss factor in structural damping is identical to the damping ratio in all analyses.", False),
    (31, "In the principle of virtual work for static equilibrium, the effects of internal forces and ideal constraint forces are eliminated.", True),
    (32, "Virtual displacement is an actual displacement that the system undergoes during motion.", False),
    (33, "D'Alembert's principle can be interpreted as the principle of virtual work applied to dynamic equilibrium.", True),
    (34, "Hamilton's principle can be written as delta integral of L dt = 0 for conservative holonomic systems.", True),
    (35, "Lagrange's equation can be derived only from Newton's second law and cannot be derived from Hamilton's principle.", False),
    (36, "The number of degrees of freedom is the minimum number of independent coordinates required to describe the configuration completely.", True),
    (37, "Generalized coordinates in this course are required to be independent.", True),
    (38, "If the number of geometric constraints increases while the number of particles is fixed, the number of DOF increases.", False),
    (39, "An equilibrium point in state space means that the state is constant in time.", True),
    (40, "A trivial equilibrium point is the only possible equilibrium point for every nonlinear vibration system.", False),
    (41, "If there exists a positive definite Lyapunov function whose time derivative is negative definite, then the equilibrium point is asymptotically stable.", True),
    (42, "If there exists a positive definite Lyapunov function whose time derivative is negative semidefinite, then the equilibrium point is necessarily asymptotically stable.", False),
    (43, "If the total energy is positive definite and its time derivative is negative definite due to pervasive damping, the equilibrium can be asymptotically stable.", True),
    (44, "Stable and asymptotically stable mean exactly the same thing.", False),
    (45, "A linearized system with a negative effective stiffness can show divergence and be unstable.", True),
    (46, "For a real symmetric positive definite matrix, all eigenvalues are real and positive.", True),
    (47, "For a real symmetric matrix, eigenvectors corresponding to distinct eigenvalues can be chosen orthogonal.", True),
    (48, "In modal analysis of an undamped conservative MDOF system, the modal matrix can be normalized so that U^T M U = I.", True),
    (49, "In modal analysis of an undamped conservative MDOF system, the transformed stiffness matrix is generally dense and non-diagonal.", False),
    (50, "Mode shapes can be multiplied by arbitrary nonzero constants without changing the eigenvalue.", True),
    (51, "If the damping matrix is proportional to M and K, the modal equations become decoupled in modal coordinates.", True),
    (52, "If the damping matrix is non-proportional, the real undamped modes always decouple the equations exactly.", False),
    (53, "For a stable linear state equation, all eigenvalues must have negative real parts.", False),
    (54, "For an asymptotically stable linear system, every eigenvalue has nonpositive real part and at least one eigenvalue may have positive real part.", False),
    (55, "Complex conjugate eigenvalues with negative real parts correspond to decaying oscillatory motion.", True),
    (56, "In a gyroscopic conservative system, the gyroscopic matrix is skew-symmetric.", True),
    (57, "A gyroscopic conservative system must always be unstable because its system matrix is asymmetric.", False),
    (58, "In a gyroscopic system, forward whirl mode and backward whirl mode can appear.", True),
    (59, "For a nonsymmetric eigenvalue problem, right eigenvectors alone generally do not satisfy ordinary orthogonality, and left eigenvectors are needed for bi-orthogonality.", True),
    (60, "A gyroscopic asymmetric system cannot be diagonalized in modal analysis using right and left eigenvectors.", False),
]

correct = 0
wrong = []
for num, stmt, expected in questions:
    result = solver.solve({"statement": stmt})
    verdict = result.final_answer.get("verdict", "N/A")
    expected_str = "TRUE" if expected else "FALSE"
    if verdict == expected_str:
        correct += 1
    else:
        rule_id = result.final_answer.get("rule_id", "?")
        confidence = result.final_answer.get("confidence", "?")
        wrong.append((num, stmt[:80], expected_str, verdict, rule_id, confidence))

print(f"Result: {correct}/60 correct ({correct/60*100:.1f}%)")
print(f"Wrong: {len(wrong)} questions")
print()
if wrong:
    for num, stmt, exp, got, rule, conf in wrong:
        print(f"  Q{num:2d} X  expected={exp:5s} got={got:5s}  rule={rule}  conf={conf}")
        print(f"       {stmt}...")
        print()
