# ME551 Additional True/False Problem Set (100 Questions)
이 문서는 업로드된 ME551 강의노트와 기출문제의 스타일을 반영해 만든 추가 T/F 문제집입니다.
구성은 선형시스템, 상태방정식, 동역학 원리, SDOF 응답, 평형점·안정성, 모달해석, 자이로·비보존계까지 고르게 분포시켰습니다.
각 문항은 True/False와 한 줄 이유를 함께 넣어 복습용으로 바로 사용할 수 있게 정리했습니다.

## 1. Linear System Theory, Convolution, FRF

1. A linear differential operator D satisfies D[αx1+βx2]=αD[x1]+βD[x2] for arbitrary constants α, β.
   Answer: T
   Reason: This is exactly the definition of linearity, combining homogeneity and additivity.

2. If a system is time-invariant, then a delayed input always produces a delayed output of the same shape.
   Answer: T
   Reason: Time invariance means a time shift at the input causes the same time shift at the output.

3. Every time-invariant system is linear.
   Answer: F
   Reason: Time invariance and linearity are different properties; a system may satisfy one without satisfying the other.

4. For an LTI system, the total response can be decomposed into a homogeneous part and a particular part.
   Answer: T
   Reason: This is the standard decomposition into free response and forced response.

5. An impulsive excitation in a second-order system is equivalent to an initial displacement jump, not an initial velocity jump.
   Answer: F
   Reason: For a second-order system, a force impulse is equivalent to an initial velocity change.

6. For a first-order system cx_dot + kx = δ(t), the impulse response is proportional to e^(-(k/c)t) for t>0.
   Answer: T
   Reason: The first-order impulse response decays exponentially with time constant c/k.

7. For an LTI system, the convolution integral is used only when the input is harmonic.
   Answer: F
   Reason: Convolution applies to arbitrary inputs, not only harmonic ones.

8. In the frequency domain, convolution in time becomes multiplication.
   Answer: T
   Reason: This is a core transform property used in FRF-based response analysis.

9. If X(jω)=G(jω)F(jω), then G(jω) is called the compliance FRF when output is displacement and input is force.
   Answer: T
   Reason: Displacement over force is the compliance or dynamic flexibility FRF.

10. The unit step function is the derivative of the Dirac delta function.
   Answer: F
   Reason: The Dirac delta is the derivative of the unit step, not the other way around.

11. The unit ramp function is the time integral of the unit step function.
   Answer: T
   Reason: Since r(t)=∫u(τ)dτ, the ramp is the integral of the step.

12. Laplace transform methods require zero initial conditions in every vibration problem.
   Answer: F
   Reason: Zero initial conditions simplify the transfer-function form, but Laplace transforms can also include nonzero initial conditions.

13. A transfer function G(s)=X(s)/F(s) is defined most cleanly under zero initial conditions.
   Answer: T
   Reason: That is the standard definition used in linear system theory.

14. The poles of G(s) are the roots of the numerator polynomial.
   Answer: F
   Reason: Poles are roots of the denominator; numerator roots are zeros.

15. If all poles of a scalar continuous-time LTI system have strictly negative real parts, the zero-input response decays to zero.
   Answer: T
   Reason: Negative real parts imply asymptotic decay of all modal terms.

16. The impulse response of an LTI system fully characterizes its input-output behavior.
   Answer: T
   Reason: Any response can be obtained from the impulse response by convolution.

17. For a causal LTI system, the impulse response must be nonzero for negative time.
   Answer: F
   Reason: Causality requires the impulse response to be zero for negative time.

18. A harmonic input to a linear system produces a steady-state output at the same frequency.
   Answer: T
   Reason: Linearity preserves frequency for sinusoidal steady-state response.

19. A harmonic input to a nonlinear system must also produce only that same frequency.
   Answer: F
   Reason: Nonlinearities generally generate harmonics and combination frequencies.

20. In Parseval-type relations, response energy can be represented equivalently in the time domain or frequency domain.
   Answer: T
   Reason: The theorem expresses the same energy measure in both domains.

## 2. State Equation and Transition Matrix

21. A p-th order linear system with n outputs can be rewritten as a first-order state equation of dimension pn.
   Answer: T
   Reason: State augmentation stacks variables and derivatives into a first-order form.

22. A second-order N-DOF vibration system becomes a first-order state system of dimension N after augmentation.
   Answer: F
   Reason: Its state dimension becomes 2N.

23. For x_dot = Ax, the transition matrix is e^{At}.
   Answer: T
   Reason: This is the exact homogeneous solution operator.

24. The transition matrix depends on the initial state x(0).
   Answer: F
   Reason: e^{At} depends only on A and time, not on a specific initial state.

25. For x_dot = Ax + Bf(t), the forced response under x(0)=0 is x(t)=∫_0^t e^{A(t-τ)}Bf(τ)dτ.
   Answer: T
   Reason: This is the standard convolution form in state space.

26. In a state-space model, the number of states always equals the number of physical coordinates.
   Answer: F
   Reason: States may include velocities or higher derivatives, so the dimension is often larger.

27. Different state choices can represent the same physical system.
   Answer: T
   Reason: State-space representation is not unique.

28. For an LTI state matrix A, the eigenvalues of A control stability of the linearized motion.
   Answer: T
   Reason: The sign of the real parts governs decay, neutrality, or growth.

29. If one eigenvalue of A has positive real part, the equilibrium is asymptotically stable.
   Answer: F
   Reason: A positive real part makes the equilibrium unstable.

30. Purely imaginary eigenvalues in a linear system always imply asymptotic stability.
   Answer: F
   Reason: They indicate neutral or marginal behavior, not asymptotic decay.

31. A real second-order vibration problem written in state form may have a nonsymmetric A even if M, C, K are symmetric.
   Answer: T
   Reason: The augmented state matrix is generally nonsymmetric.

32. The matrix exponential e^{At} can be obtained efficiently from eigen-decomposition when A is diagonalizable.
   Answer: T
   Reason: Diagonalization converts the exponential into exponentials of eigenvalues.

33. For any square matrix, similarity transformation changes the eigenvalues.
   Answer: F
   Reason: Similarity transformations preserve eigenvalues.

34. Orthogonal transformation is a special case of similarity transformation.
   Answer: T
   Reason: It is similarity with P^{-1}=P^T.

35. If two matrices are similar, they must have identical eigenvectors in the original coordinates.
   Answer: F
   Reason: They have the same eigenvalues, but eigenvectors change with coordinates.

## 3. Principles of Dynamics and Analytical Mechanics

36. The number of degrees of freedom is the minimum number of independent generalized coordinates needed to describe the configuration.
   Answer: T
   Reason: That is the standard definition of DOF.

37. Generalized coordinates must always coincide with Cartesian coordinates.
   Answer: F
   Reason: They may be angular, abstract, or otherwise non-Cartesian.

38. A virtual displacement is an actual infinitesimal motion that necessarily satisfies Newton's second law.
   Answer: F
   Reason: It is a compatible imagined variation, not the actual motion.

39. In the principle of virtual work for ideal constraints, constraint forces do no virtual work.
   Answer: T
   Reason: They are orthogonal to admissible virtual displacements.

40. Internal action-reaction forces can be eliminated in the virtual-work summation for a system of particles.
   Answer: T
   Reason: They cancel in pairs under ideal assumptions.

41. D'Alembert's principle can be interpreted as virtual work applied to dynamic equilibrium.
   Answer: T
   Reason: It introduces inertia forces so the dynamic problem resembles static equilibrium.

42. The generalized inertial force in D'Alembert's principle is treated with the same sign as the applied generalized force.
   Answer: F
   Reason: It enters with opposite sign when moved into the virtual-work expression.

43. Hamilton's principle uses path variations that vanish at the two endpoint times.
   Answer: T
   Reason: The endpoint variations are prescribed to be zero.

44. For conservative holonomic systems, Hamilton's principle is δ∫L dt = 0.
   Answer: T
   Reason: This is the standard variational statement.

45. Lagrange's equation is d/dt(∂L/∂q_dot)-∂L/∂q = Q_nc for generalized nonconservative force Q_nc.
   Answer: T
   Reason: This is the basic form used throughout analytical dynamics.

46. A system with a viscous damper can always be represented using potential energy only.
   Answer: F
   Reason: Damping is dissipative and is not captured by ordinary potential energy.

47. Kinetic energy is always quadratic in generalized speeds for standard rigid or particle systems.
   Answer: T
   Reason: It has the general quadratic form in q_dot.

48. For a natural system in the lecture-note sense, the coordinate frame is not rotating relative to the inertial frame.
   Answer: T
   Reason: Then gyroscopic and centrifugal-type terms do not appear from coordinate rotation.

49. In generalized coordinates, independence of coordinates is essential in deriving separate equations from the virtual-work statement.
   Answer: T
   Reason: The arbitrariness of each δq_k is what yields the individual equations.

50. A non-holonomic constraint is one that depends only on positions and time.
   Answer: F
   Reason: That description is holonomic; non-holonomic constraints may involve velocities.

## 4. Single-DOF Systems and Harmonic Response

51. For the free response of an undamped SDOF system, the motion is purely harmonic with frequency ω_n = √(k/m).
   Answer: T
   Reason: The solution is sinusoidal at the natural frequency.

52. For an underdamped SDOF system, the damped natural frequency ω_d is greater than ω_n.
   Answer: F
   Reason: ω_d = ω_n√(1-ζ^2), so it is smaller when 0<ζ<1.

53. Critical damping corresponds to ζ = 1.
   Answer: T
   Reason: That is the defining value separating oscillatory and nonoscillatory decay.

54. If ζ>1, the free response is oscillatory but slowly decaying.
   Answer: F
   Reason: That is overdamping; the response is nonoscillatory.

55. The quality factor of a lightly damped second-order system is approximately 1/(2ζ).
   Answer: T
   Reason: This is a standard small-damping approximation.

56. At resonance in an undamped harmonically forced SDOF system, the steady-state amplitude becomes unbounded.
   Answer: T
   Reason: The ideal linear undamped model predicts infinite amplitude at exact resonance.

57. Increasing viscous damping generally lowers the resonance peak in the compliance FRF.
   Answer: T
   Reason: Damping reduces dynamic amplification.

58. For a first-order system, the break frequency is equal to the inverse time constant.
   Answer: T
   Reason: ω_b = 1/τ.

59. In a first-order compliance FRF, the phase goes from 0 to -π as frequency increases.
   Answer: F
   Reason: For first order it goes from 0 toward -π/2.

60. In a second-order compliance FRF, the low-frequency limit equals the static compliance 1/k.
   Answer: T
   Reason: At ω→0 the inertial and damping terms vanish.

61. For an underdamped second-order system, the half-power points are used to estimate damping.
   Answer: T
   Reason: They are the basis of the bandwidth method.

62. Structural damping is modeled most naturally as an equivalent complex stiffness in harmonic analysis.
   Answer: T
   Reason: That is the common frequency-domain representation.

63. Structural damping is generally as reliable for transient time-domain analysis as for harmonic analysis.
   Answer: F
   Reason: The complex-stiffness model is mainly appropriate for harmonic response.

64. Equivalent viscous damping is introduced by matching energy dissipation per cycle.
   Answer: T
   Reason: That is the standard equivalence criterion.

65. In viscous damping, energy dissipation per cycle is independent of frequency for a fixed amplitude.
   Answer: F
   Reason: For viscous damping, cycle loss depends on frequency.

66. A force impulse on a first-order system is equivalent to an initial displacement change.
   Answer: T
   Reason: The impulse causes a state jump in the first-order variable.

67. A force impulse on a second-order system is equivalent to an initial velocity change.
   Answer: T
   Reason: The momentum jump appears as a change in x_dot.

68. The unit impulse response of a stable first-order system is nondecaying periodic motion.
   Answer: F
   Reason: It is exponentially decaying, not periodic.

69. For zero initial conditions, multiplying G(s) by 1/s corresponds to the response to unit-step input.
   Answer: T
   Reason: Since L{u(t)}=1/s.

70. A periodic excitation can be expanded in a Fourier series and the response found by superposing responses to each harmonic.
   Answer: T
   Reason: This is the standard periodic forcing approach in linear systems.

## 5. Equilibrium, Stability, and Modal Analysis

71. An equilibrium point in state space is a point at which the state remains constant in time.
   Answer: T
   Reason: That is the definition of equilibrium.

72. Every equilibrium point is the origin.
   Answer: F
   Reason: The origin is only the trivial equilibrium; nontrivial equilibria may also exist.

73. In configuration space alone, different dynamical paths can intersect even when the actual state trajectories are distinct.
   Answer: T
   Reason: Because configuration space omits velocity information.

74. State-space trajectories for a deterministic system with different initial conditions cannot intersect in finite time.
   Answer: T
   Reason: Otherwise uniqueness of solutions would be violated.

75. For determining equilibrium in Lagrangian systems, one sets q_dot=q_ddot=0 in the equations of motion.
   Answer: T
   Reason: Equilibrium requires constant generalized coordinates.

76. Potential energy must have a stationary value at an equilibrium of a natural conservative system.
   Answer: T
   Reason: Equilibrium follows from vanishing gradient of the potential.

77. If the linearized stiffness about an equilibrium is negative in an SDOF problem, the equilibrium is stable.
   Answer: F
   Reason: Negative effective stiffness leads to exponential growth and instability.

78. A Lyapunov function must be unique if it exists.
   Answer: F
   Reason: Many different Lyapunov functions may work for the same equilibrium.

79. If a positive definite Lyapunov function has negative definite time derivative, the equilibrium is asymptotically stable.
   Answer: T
   Reason: This is the standard Lyapunov theorem.

80. If both the Lyapunov function and its time derivative are sign-variable, no definite stability conclusion follows from that test alone.
   Answer: T
   Reason: The direct method becomes inconclusive in that case.

81. For a real symmetric eigenvalue problem Ku=λMu with M positive definite, the eigenvalues λ are real.
   Answer: T
   Reason: The transformed standard problem is symmetric.

82. For a conservative natural system, mode shapes can be mass-normalized so that U^T M U = I.
   Answer: T
   Reason: Mass normalization is standard in modal analysis.

83. For a conservative natural system, modal coordinates decouple the undamped equations of motion.
   Answer: T
   Reason: That is the purpose of modal expansion.

84. If damping is proportional, the same undamped mode shapes can be used to decouple the damped equations.
   Answer: T
   Reason: Classical damping preserves modal decoupling.

85. Any nonproportional damping matrix can still be decoupled by the real undamped mode matrix alone.
   Answer: F
   Reason: General nonclassical damping usually does not decouple that way.

86. In modal analysis, multiplying a mode shape by a nonzero scalar changes the natural frequency.
   Answer: F
   Reason: Eigenvectors are scale-indeterminate; eigenvalues do not change.

87. Repeated eigenvalues in symmetric problems may require special care in choosing an orthogonal modal basis.
   Answer: T
   Reason: Within the repeated eigenspace, basis vectors are not unique.

88. For a symmetric matrix, orthogonality of distinct eigenvectors follows automatically when eigenvalues are distinct.
   Answer: T
   Reason: This is a fundamental result of the symmetric EVP.

89. The modal matrix of a conservative system diagonalizes both the mass and stiffness matrices after proper normalization.
   Answer: T
   Reason: It yields U^T M U = I and U^T K U = Λ.

90. In a lightly damped modal equation, each modal coordinate behaves like an independent damped SDOF oscillator.
   Answer: T
   Reason: That is exactly the decoupled modal interpretation.

## 6. Gyroscopic and Nonconservative Systems

91. A gyroscopic matrix G in Mq_ddot + Gq_dot + Kq = 0 is symmetric.
   Answer: F
   Reason: The gyroscopic matrix is skew-symmetric.

92. Gyroscopic systems can have complex mode shapes even when the physical coordinates are real.
   Answer: T
   Reason: The physical solution is real but individual mode vectors can be complex.

93. A conservative gyroscopic system must be unstable because its modes are complex.
   Answer: F
   Reason: Complex mode shapes do not by themselves imply instability.

94. Forward whirl and backward whirl can appear in rotating systems with gyroscopic effects.
   Answer: T
   Reason: These are characteristic rotating-mode phenomena.

95. When rotation speed Ω goes to zero in a gyroscopic system, the x and y motions may reduce to ordinary decoupled natural motions.
   Answer: T
   Reason: The gyroscopic coupling disappears when Ω=0.

96. In a nonsymmetric eigenvalue problem, right eigenvectors alone generally satisfy the same orthogonality relation as in the symmetric case.
   Answer: F
   Reason: Standard orthogonality is lost; left-right biorthogonality is used instead.

97. The eigenvalues of A and A^T are the same.
   Answer: T
   Reason: Transpose preserves the characteristic polynomial.

98. For a nonsymmetric matrix, left eigenvectors are obtained from the adjoint problem A^T y = λ y.
   Answer: T
   Reason: These provide the dual basis used in biorthogonality.

99. Biorthogonality means y_i^T x_j = 0 for i≠j after suitable normalization.
   Answer: T
   Reason: This is the key replacement for orthogonality in nonsymmetric systems.

100. If a linearized system has an eigenvalue with positive real part and nonzero imaginary part, that always indicates a decaying oscillation.
   Answer: F
   Reason: Positive real part means growing oscillation, i.e., flutter-type instability.

