"""Verify ConceptDBSolver against 100 extra T/F questions."""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
from me551_solver.core.concept_db import ConceptDBSolver

solver = ConceptDBSolver()

questions = [
    (1, "A linear differential operator D satisfies D[alpha*x1+beta*x2]=alpha*D[x1]+beta*D[x2] for arbitrary constants alpha, beta.", True),
    (2, "If a system is time-invariant, then a delayed input always produces a delayed output of the same shape.", True),
    (3, "Every time-invariant system is linear.", False),
    (4, "For an LTI system, the total response can be decomposed into a homogeneous part and a particular part.", True),
    (5, "An impulsive excitation in a second-order system is equivalent to an initial displacement jump, not an initial velocity jump.", False),
    (6, "For a first-order system cx_dot + kx = delta(t), the impulse response is proportional to e^(-(k/c)t) for t>0.", True),
    (7, "For an LTI system, the convolution integral is used only when the input is harmonic.", False),
    (8, "In the frequency domain, convolution in time becomes multiplication.", True),
    (9, "If X(jw)=G(jw)F(jw), then G(jw) is called the compliance FRF when output is displacement and input is force.", True),
    (10, "The unit step function is the derivative of the Dirac delta function.", False),
    (11, "The unit ramp function is the time integral of the unit step function.", True),
    (12, "Laplace transform methods require zero initial conditions in every vibration problem.", False),
    (13, "A transfer function G(s)=X(s)/F(s) is defined most cleanly under zero initial conditions.", True),
    (14, "The poles of G(s) are the roots of the numerator polynomial.", False),
    (15, "If all poles of a scalar continuous-time LTI system have strictly negative real parts, the zero-input response decays to zero.", True),
    (16, "The impulse response of an LTI system fully characterizes its input-output behavior.", True),
    (17, "For a causal LTI system, the impulse response must be nonzero for negative time.", False),
    (18, "A harmonic input to a linear system produces a steady-state output at the same frequency.", True),
    (19, "A harmonic input to a nonlinear system must also produce only that same frequency.", False),
    (20, "In Parseval-type relations, response energy can be represented equivalently in the time domain or frequency domain.", True),
    (21, "A p-th order linear system with n outputs can be rewritten as a first-order state equation of dimension pn.", True),
    (22, "A second-order N-DOF vibration system becomes a first-order state system of dimension N after augmentation.", False),
    (23, "For x_dot = Ax, the transition matrix is e^{At}.", True),
    (24, "The transition matrix depends on the initial state x(0).", False),
    (25, "For x_dot = Ax + Bf(t), the forced response under x(0)=0 is x(t)=integral of e^{A(t-tau)}Bf(tau)dtau.", True),
    (26, "In a state-space model, the number of states always equals the number of physical coordinates.", False),
    (27, "Different state choices can represent the same physical system.", True),
    (28, "For an LTI state matrix A, the eigenvalues of A control stability of the linearized motion.", True),
    (29, "If one eigenvalue of A has positive real part, the equilibrium is asymptotically stable.", False),
    (30, "Purely imaginary eigenvalues in a linear system always imply asymptotic stability.", False),
    (31, "A real second-order vibration problem written in state form may have a nonsymmetric A even if M, C, K are symmetric.", True),
    (32, "The matrix exponential e^{At} can be obtained efficiently from eigen-decomposition when A is diagonalizable.", True),
    (33, "For any square matrix, similarity transformation changes the eigenvalues.", False),
    (34, "Orthogonal transformation is a special case of similarity transformation.", True),
    (35, "If two matrices are similar, they must have identical eigenvectors in the original coordinates.", False),
    (36, "The number of degrees of freedom is the minimum number of independent generalized coordinates needed to describe the configuration.", True),
    (37, "Generalized coordinates must always coincide with Cartesian coordinates.", False),
    (38, "A virtual displacement is an actual infinitesimal motion that necessarily satisfies Newton's second law.", False),
    (39, "In the principle of virtual work for ideal constraints, constraint forces do no virtual work.", True),
    (40, "Internal action-reaction forces can be eliminated in the virtual-work summation for a system of particles.", True),
    (41, "D'Alembert's principle can be interpreted as virtual work applied to dynamic equilibrium.", True),
    (42, "The generalized inertial force in D'Alembert's principle is treated with the same sign as the applied generalized force.", False),
    (43, "Hamilton's principle uses path variations that vanish at the two endpoint times.", True),
    (44, "For conservative holonomic systems, Hamilton's principle is delta integral of L dt = 0.", True),
    (45, "Lagrange's equation is d/dt(dL/dq_dot)-dL/dq = Q_nc for generalized nonconservative force Q_nc.", True),
    (46, "A system with a viscous damper can always be represented using potential energy only.", False),
    (47, "Kinetic energy is always quadratic in generalized speeds for standard rigid or particle systems.", True),
    (48, "For a natural system in the lecture-note sense, the coordinate frame is not rotating relative to the inertial frame.", True),
    (49, "In generalized coordinates, independence of coordinates is essential in deriving separate equations from the virtual-work statement.", True),
    (50, "A non-holonomic constraint is one that depends only on positions and time.", False),
    (51, "For the free response of an undamped SDOF system, the motion is purely harmonic with frequency omega_n = sqrt(k/m).", True),
    (52, "For an underdamped SDOF system, the damped natural frequency omega_d is greater than omega_n.", False),
    (53, "Critical damping corresponds to zeta = 1.", True),
    (54, "If zeta>1, the free response is oscillatory but slowly decaying.", False),
    (55, "The quality factor of a lightly damped second-order system is approximately 1/(2*zeta).", True),
    (56, "At resonance in an undamped harmonically forced SDOF system, the steady-state amplitude becomes unbounded.", True),
    (57, "Increasing viscous damping generally lowers the resonance peak in the compliance FRF.", True),
    (58, "For a first-order system, the break frequency is equal to the inverse time constant.", True),
    (59, "In a first-order compliance FRF, the phase goes from 0 to -pi as frequency increases.", False),
    (60, "In a second-order compliance FRF, the low-frequency limit equals the static compliance 1/k.", True),
    (61, "For an underdamped second-order system, the half-power points are used to estimate damping.", True),
    (62, "Structural damping is modeled most naturally as an equivalent complex stiffness in harmonic analysis.", True),
    (63, "Structural damping is generally as reliable for transient time-domain analysis as for harmonic analysis.", False),
    (64, "Equivalent viscous damping is introduced by matching energy dissipation per cycle.", True),
    (65, "In viscous damping, energy dissipation per cycle is independent of frequency for a fixed amplitude.", False),
    (66, "A force impulse on a first-order system is equivalent to an initial displacement change.", True),
    (67, "A force impulse on a second-order system is equivalent to an initial velocity change.", True),
    (68, "The unit impulse response of a stable first-order system is nondecaying periodic motion.", False),
    (69, "For zero initial conditions, multiplying G(s) by 1/s corresponds to the response to unit-step input.", True),
    (70, "A periodic excitation can be expanded in a Fourier series and the response found by superposing responses to each harmonic.", True),
    (71, "An equilibrium point in state space is a point at which the state remains constant in time.", True),
    (72, "Every equilibrium point is the origin.", False),
    (73, "In configuration space alone, different dynamical paths can intersect even when the actual state trajectories are distinct.", True),
    (74, "State-space trajectories for a deterministic system with different initial conditions cannot intersect in finite time.", True),
    (75, "For determining equilibrium in Lagrangian systems, one sets q_dot=q_ddot=0 in the equations of motion.", True),
    (76, "Potential energy must have a stationary value at an equilibrium of a natural conservative system.", True),
    (77, "If the linearized stiffness about an equilibrium is negative in an SDOF problem, the equilibrium is stable.", False),
    (78, "A Lyapunov function must be unique if it exists.", False),
    (79, "If a positive definite Lyapunov function has negative definite time derivative, the equilibrium is asymptotically stable.", True),
    (80, "If both the Lyapunov function and its time derivative are sign-variable, no definite stability conclusion follows from that test alone.", True),
    (81, "For a real symmetric eigenvalue problem Ku=lambda*Mu with M positive definite, the eigenvalues lambda are real.", True),
    (82, "For a conservative natural system, mode shapes can be mass-normalized so that U^T M U = I.", True),
    (83, "For a conservative natural system, modal coordinates decouple the undamped equations of motion.", True),
    (84, "If damping is proportional, the same undamped mode shapes can be used to decouple the damped equations.", True),
    (85, "Any nonproportional damping matrix can still be decoupled by the real undamped mode matrix alone.", False),
    (86, "In modal analysis, multiplying a mode shape by a nonzero scalar changes the natural frequency.", False),
    (87, "Repeated eigenvalues in symmetric problems may require special care in choosing an orthogonal modal basis.", True),
    (88, "For a symmetric matrix, orthogonality of distinct eigenvectors follows automatically when eigenvalues are distinct.", True),
    (89, "The modal matrix of a conservative system diagonalizes both the mass and stiffness matrices after proper normalization.", True),
    (90, "In a lightly damped modal equation, each modal coordinate behaves like an independent damped SDOF oscillator.", True),
    (91, "A gyroscopic matrix G in Mq_ddot + Gq_dot + Kq = 0 is symmetric.", False),
    (92, "Gyroscopic systems can have complex mode shapes even when the physical coordinates are real.", True),
    (93, "A conservative gyroscopic system must be unstable because its modes are complex.", False),
    (94, "Forward whirl and backward whirl can appear in rotating systems with gyroscopic effects.", True),
    (95, "When rotation speed Omega goes to zero in a gyroscopic system, the x and y motions may reduce to ordinary decoupled natural motions.", True),
    (96, "In a nonsymmetric eigenvalue problem, right eigenvectors alone generally satisfy the same orthogonality relation as in the symmetric case.", False),
    (97, "The eigenvalues of A and A^T are the same.", True),
    (98, "For a nonsymmetric matrix, left eigenvectors are obtained from the adjoint problem A^T y = lambda y.", True),
    (99, "Biorthogonality means y_i^T x_j = 0 for i!=j after suitable normalization.", True),
    (100, "If a linearized system has an eigenvalue with positive real part and nonzero imaginary part, that always indicates a decaying oscillation.", False),
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

print(f"Result: {correct}/100 correct ({correct}%)")
print(f"Wrong: {len(wrong)} questions")
print()
if wrong:
    for num, stmt, exp, got, rule, conf in wrong:
        print(f"  Q{num:3d} X  expected={exp:5s} got={got:5s}  rule={rule}  conf={conf}")
        print(f"       {stmt}...")
        print()
