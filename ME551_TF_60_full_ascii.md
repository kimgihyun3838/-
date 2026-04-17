# ME551 True/False Problem Set (60 Questions)
## Based on lecture notes and past midterm styles

### 사용 안내
- 각 문항에 대해 True 또는 False를 판단한 뒤, 바로 아래의 정답과 이유를 확인하십시오.
- 일부 문항은 의도적으로 기출 스타일의 함정을 포함하고 있습니다.
- 특히 다음을 자주 구분하도록 구성했습니다.
  - free vibration vs forced vibration
  - convolution in time domain vs multiplication in frequency domain
  - N차원 state equation vs 2N차원 state equation
  - symmetric EVP vs nonsymmetric EVP
  - proportional damping vs non-proportional damping
  - structural damping의 적용 범위
  - stable vs asymptotically stable

---

## 1. Linear systems, convolution, FRF

### 1
A linear system satisfies homogeneity and additivity, and therefore satisfies the superposition principle.

**Answer:** True  
**Reason:** 선형성의 정의가 바로 homogeneity와 additivity를 동시에 만족하는 것이다. 따라서 superposition principle이 성립한다.

### 2
A nonlinear system must violate both homogeneity and additivity simultaneously.

**Answer:** False  
**Reason:** 둘 다 깨질 필요는 없다. 둘 중 하나만 성립하지 않아도 nonlinear system이다.

### 3
For an arbitrary excitation of an LTI system, the response can be obtained by convolution of the input with the impulse response in the time domain.

**Answer:** True  
**Reason:** LTI 시스템에서는 \(x(t)=\int_0^t g(t-\tau)f(\tau)\,d\tau\) 형태의 convolution integral이 성립한다.

### 4
For an arbitrary excitation of an LTI system, the response in the frequency domain is obtained by convolution of the FRF and the input spectrum.

**Answer:** False  
**Reason:** frequency domain에서는 convolution이 아니라 multiplication이다. 즉 \(X(j\omega)=G(j\omega)F(j\omega)\) 이다.

### 5
The FRF method and the convolution method can be used to obtain the response due to initial conditions only.

**Answer:** False  
**Reason:** initial-condition response는 free vibration 문제이다. convolution과 FRF는 forcing term이 있을 때의 forced response를 다루는 도구다. initial conditions만 있는 경우에는 homogeneous solution을 사용해야 한다.

---

## 2. State equation and transition matrix

### 6
A second-order linear vibration system with N degrees of freedom can be rewritten as a first-order state equation of dimension 2N.

**Answer:** True  
**Reason:** 상태변수를 일반적으로 \(x=[q^T\ \dot q^T]^T\) 로 잡으면 차원이 2N이 된다.

### 7
A second-order linear MIMO system with N degrees of freedom can always be rewritten as a first-order state equation of dimension N.

**Answer:** False  
**Reason:** 일반적으로 2N차원이 필요하다. 이것이 2025 기출 T/F의 대표 함정이다.

### 8
For a p-th order linear MIMO system with n degrees of freedom, the first-order state equation generally has dimension pn.

**Answer:** True  
**Reason:** 각 DOF마다 p차 미분이 포함되므로 보통 상태변수 차원은 \(pn\) 이 된다.

### 9
The transition matrix of a first-order state equation \(\dot x=Ax\) is \(e^{At}\).

**Answer:** True  
**Reason:** homogeneous solution은 \(x(t)=e^{At}x(0)\) 이다.

### 10
The transition matrix changes the eigenvalues of the system matrix A.

**Answer:** False  
**Reason:** \(e^{At}\) 는 A로부터 만들어진 matrix exponential이며, A의 고유구조를 바탕으로 homogeneous response를 표현한다. similarity transformation이 아닌 자체 상태전이 표현이며, 원래 A의 stability 정보는 보존된다.

---

## 3. Generalized functions and impulse/step/ramp responses

### 11
The unit step function is the time integral of the Dirac delta function.

**Answer:** True  
**Reason:** \(u(t-a)=\int_{-\infty}^{t}\delta(\tau-a)\,d\tau\) 로 볼 수 있다.

### 12
The unit ramp function is the derivative of the unit step function.

**Answer:** False  
**Reason:** 반대다. ramp는 step의 적분이다. step의 미분이 impulse다.

### 13
For a first-order system, a unit impulse input is mathematically equivalent to an initial displacement of magnitude \(1/c\).

**Answer:** True  
**Reason:** 강의노트에서 1차계 \(cẋ+kx=\delta(t)\) 는 적절한 jump condition을 통해 \(x(0)=1/c\) 인 initial condition 문제와 동등하게 해석된다.

### 14
For a second-order system, a unit impulse input is mathematically equivalent to an initial velocity of magnitude \(1/m\).

**Answer:** True  
**Reason:** \(mẍ+cẋ+kx=\delta(t)\) 의 impulse는 속도에 jump를 만들며 \(ẋ(0^+)=1/m\) 과 동등하다.

### 15
For a second-order system, a unit impulse input is mathematically equivalent to an initial displacement of magnitude \(1/k\).

**Answer:** False  
**Reason:** impulse는 displacement jump보다 velocity jump를 만든다. 따라서 \(1/k\) 와 연결되는 것이 아니다.

---

## 4. First-order and second-order harmonic response

### 16
A first-order system is typically characterized by a break frequency or time constant.

**Answer:** True  
**Reason:** 1차계의 전형적 특성량은 \(\omega_b=k/c\) 와 같은 break frequency 또는 \(\tau=c/k\) 같은 time constant다.

### 17
An underdamped second-order system is typically characterized by the undamped natural frequency and the damping ratio.

**Answer:** True  
**Reason:** \(\omega_n=\sqrt{k/m}\), \(\zeta=c/(2\sqrt{mk})\) 로 정리한다.

### 18
For an underdamped second-order system, the damped natural frequency is larger than the undamped natural frequency.

**Answer:** False  
**Reason:** \(\omega_d=\omega_n\sqrt{1-\zeta^2}\) 이므로 \(\zeta<1\) 일 때 항상 \(\omega_d<\omega_n\) 이다.

### 19
For a harmonic input applied to a linear system, the steady-state output has the same frequency as the input.

**Answer:** True  
**Reason:** 선형시스템에서는 sinusoidal input에 대해 sinusoidal steady-state output이 같은 forcing frequency를 가진다.

### 20
For a harmonic input applied to a nonlinear system, the steady-state output must have only the input frequency.

**Answer:** False  
**Reason:** nonlinear system은 일반적으로 harmonic generation을 일으켜 배수 주파수나 조합 주파수가 생길 수 있다.

---

## 5. Dynamic stiffness, compliance, resonance

### 21
The compliance FRF is the reciprocal of the dynamic stiffness for an SDOF system.

**Answer:** True  
**Reason:** 강의노트의 정의에 따르면 \(G(s)=1/Z(s)\) 이다.

### 22
At zero frequency, the compliance of a second-order SDOF system is \(1/k\).

**Answer:** True  
**Reason:** \(G(j0)=1/k\) 이다. 정적 입력에 대한 정적 처짐과 같다.

### 23
At very high frequency, the compliance magnitude of a second-order SDOF system tends to zero.

**Answer:** True  
**Reason:** 고주파에서 inertia term이 지배적이므로 \(G(j\omega)\sim 1/(m\omega^2)\to 0\).

### 24
Increasing damping ratio generally lowers the resonance peak of the compliance FRF.

**Answer:** True  
**Reason:** damping이 커질수록 peak amplification이 감소한다.

### 25
At resonance of an undamped second-order system under harmonic force, the steady-state amplitude remains finite.

**Answer:** False  
**Reason:** 이상적인 무감쇠 계에서 exact resonance이면 steady-state amplitude가 무한대로 발산한다.

---

## 6. Structural damping

### 26
Structural damping is commonly associated with hysteretic energy loss of elastic materials under cyclic stress.

**Answer:** True  
**Reason:** 강의노트에서 structural damping의 물리적 배경을 hysteresis로 설명한다.

### 27
Structural damping is usually represented by an equivalent viscous damping only for harmonic analysis.

**Answer:** True  
**Reason:** harmonic loading과 harmonic response를 가정할 때 equivalent viscous damping으로 해석하는 것이 일반적이다.

### 28
Structural damping can be used in transient analysis with the same validity as in harmonic analysis.

**Answer:** False  
**Reason:** 기출 T/F 핵심 포인트다. structural damping의 복소강성 표현은 harmonic setting에서 주로 사용되며 transient analysis에 그대로 쓰면 인과성 등 문제가 생긴다.

### 29
For structural damping, the energy dissipation per cycle is proportional to the square of the response amplitude.

**Answer:** True  
**Reason:** 강의노트에서 \( \Delta E_{cycle}\propto x_0^2 \) 형태를 사용한다.

### 30
The loss factor in structural damping is identical to the damping ratio in all analyses.

**Answer:** False  
**Reason:** loss factor \(\eta\) 와 viscous damping ratio \(\zeta\) 는 동일한 개념이 아니다. harmonic equivalence에서는 \(\eta=2\zeta\) 같은 관계를 쓸 수 있지만 항상 동일한 것은 아니다.

---

## 7. Virtual work, D’Alembert, Hamilton, Lagrange

### 31
In the principle of virtual work for static equilibrium, the effects of internal forces and ideal constraint forces are eliminated.

**Answer:** True  
**Reason:** internal forces는 action-reaction으로 상쇄되고, ideal constraint force는 admissible virtual displacement와 수직이라 virtual work를 하지 않는다.

### 32
Virtual displacement is an actual displacement that the system undergoes during motion.

**Answer:** False  
**Reason:** virtual displacement는 실제 시간발전으로 생긴 displacement가 아니라, constraints와 양립하는 imagined infinitesimal variation이다.

### 33
D’Alembert’s principle can be interpreted as the principle of virtual work applied to dynamic equilibrium.

**Answer:** True  
**Reason:** inertia force를 추가하여 dynamic problem을 equilibrium form으로 바꾸는 것이 핵심이다.

### 34
Hamilton’s principle can be written as \(\delta\int_{t_1}^{t_2}L\,dt=0\) for conservative holonomic systems.

**Answer:** True  
**Reason:** 비보존력이 없고 holonomic이면 extended term 없이 그 형태로 쓴다.

### 35
Lagrange’s equation can be derived only from Newton’s second law and cannot be derived from Hamilton’s principle.

**Answer:** False  
**Reason:** 강의노트에서는 Hamilton’s principle로부터 유도하는 흐름을 명확히 보여준다. D’Alembert principle로부터도 유도할 수 있다.

---

## 8. Generalized coordinates, DOF, equilibrium

### 36
The number of degrees of freedom is the minimum number of independent coordinates required to describe the configuration completely.

**Answer:** True  
**Reason:** DOF의 정의 그 자체다.

### 37
Generalized coordinates in this course are required to be independent.

**Answer:** True  
**Reason:** 강의노트에서 generalized coordinates는 independent of each other 라고 명시한다.

### 38
If the number of geometric constraints increases while the number of particles is fixed, the number of DOF increases.

**Answer:** False  
**Reason:** 제약이 늘면 독립좌표 수는 감소한다.

### 39
An equilibrium point in state space means that the state is constant in time.

**Answer:** True  
**Reason:** equilibrium state는 \(\dot x=0\) 을 만족하는 constant state다. 반드시 모든 generalized coordinates가 zero일 필요는 없다.

### 40
A trivial equilibrium point is the only possible equilibrium point for every nonlinear vibration system.

**Answer:** False  
**Reason:** rotating hoop, rotating disk 문제처럼 parameter에 따라 여러 nontrivial equilibrium points가 생길 수 있다.

---

## 9. Stability and Lyapunov method

### 41
If there exists a positive definite Lyapunov function whose time derivative is negative definite, then the equilibrium point is asymptotically stable.

**Answer:** True  
**Reason:** Lyapunov direct method의 대표 정리다.

### 42
If there exists a positive definite Lyapunov function whose time derivative is negative semidefinite, then the equilibrium point is necessarily asymptotically stable.

**Answer:** False  
**Reason:** negative semidefinite면 일반적으로 stable까지만 바로 말할 수 있다. asymptotically stable을 보장하려면 추가 조건이 필요하다.

### 43
If the total energy is positive definite and its time derivative is negative definite due to pervasive damping, the equilibrium can be asymptotically stable.

**Answer:** True  
**Reason:** 강의노트에서 energy-like Lyapunov function과 damping을 이용해 asymptotic stability를 설명한다.

### 44
Stable and asymptotically stable mean exactly the same thing.

**Answer:** False  
**Reason:** stable은 근처에 머무는 것이고, asymptotically stable은 시간이 갈수록 equilibrium으로 수렴하는 것이다.

### 45
A linearized system with a negative effective stiffness can show divergence and be unstable.

**Answer:** True  
**Reason:** linearized equation의 restoring term 부호가 잘못되면 exponential divergence가 발생할 수 있다. rotating triangle 예제가 대표적이다.

---

## 10. Symmetric eigenvalue problem and modal analysis

### 46
For a real symmetric positive definite matrix, all eigenvalues are real and positive.

**Answer:** True  
**Reason:** lecture note와 2025 기출 T/F의 핵심 사실이다.

### 47
For a real symmetric matrix, eigenvectors corresponding to distinct eigenvalues can be chosen orthogonal.

**Answer:** True  
**Reason:** symmetric EVP의 fundamental property이다.

### 48
In modal analysis of an undamped conservative MDOF system, the modal matrix can be normalized so that \(U^TMU=I\).

**Answer:** True  
**Reason:** mass-normalization을 하면 그렇게 둘 수 있다.

### 49
In modal analysis of an undamped conservative MDOF system, the transformed stiffness matrix is generally dense and non-diagonal.

**Answer:** False  
**Reason:** 적절한 modal matrix를 쓰면 \(U^TKU\) 는 diagonal이 된다.

### 50
Mode shapes can be multiplied by arbitrary nonzero constants without changing the eigenvalue.

**Answer:** True  
**Reason:** eigenvector는 scale-invariant하다.

---

## 11. Proportional damping, non-proportional damping, state-space eigenvalues

### 51
If the damping matrix is proportional to M and K, the modal equations become decoupled in modal coordinates.

**Answer:** True  
**Reason:** classical damping 또는 proportional damping의 핵심 장점이다.

### 52
If the damping matrix is non-proportional, the real undamped modes always decouple the equations exactly.

**Answer:** False  
**Reason:** non-proportional damping이면 일반적인 real modal transformation만으로 완전 decoupling이 되지 않는다.

### 53
For a stable linear state equation, all eigenvalues must have negative real parts.

**Answer:** False  
**Reason:** asymptotically stable이라면 negative real parts가 필요하다. 단순 stable은 purely imaginary eigenvalues를 가질 수 있다.

### 54
For an asymptotically stable linear system, every eigenvalue has nonpositive real part and at least one eigenvalue may have positive real part.

**Answer:** False  
**Reason:** positive real part가 하나라도 있으면 unstable이다.

### 55
Complex conjugate eigenvalues with negative real parts correspond to decaying oscillatory motion.

**Answer:** True  
**Reason:** \(e^{(\sigma\pm j\omega)t}\) 에서 \(\sigma<0\) 이면 진동하면서 감쇠한다.

---

## 12. Gyroscopic systems and nonsymmetric eigenvalue problems

### 56
In a gyroscopic conservative system, the gyroscopic matrix is skew-symmetric.

**Answer:** True  
**Reason:** 강의노트에서 \(G=-G^T\) 인 skew-symmetric matrix로 정의한다.

### 57
A gyroscopic conservative system must always be unstable because its system matrix is asymmetric.

**Answer:** False  
**Reason:** asymmetry 자체가 곧 instability를 뜻하지 않는다. gyroscopic conservative system은 pure imaginary eigenvalues를 가질 수 있고 stable whirl motion이 가능하다.

### 58
In a gyroscopic system, forward whirl mode and backward whirl mode can appear.

**Answer:** True  
**Reason:** lecture 12에서 rotating disk 예제를 통해 forward whirl과 backward whirl을 구분한다.

### 59
For a nonsymmetric eigenvalue problem, right eigenvectors alone generally do not satisfy ordinary orthogonality, and left eigenvectors are needed for bi-orthogonality.

**Answer:** True  
**Reason:** lecture 13의 핵심 내용이다. unsymmetric A에 대해서는 right/left eigenvectors의 bi-orthogonality를 사용한다.

### 60
A gyroscopic asymmetric system cannot be diagonalized in modal analysis using right and left eigenvectors.

**Answer:** False  
**Reason:** 2025 midterm solution의 핵심 포인트다. ordinary orthogonality는 안 되지만 right and left eigenvectors를 이용한 bi-orthogonal modal analysis는 가능하다.

---

## Quick review checklist

다음 문항들을 특별히 반복 확인하십시오.

- 4번, 5번: convolution과 multiplication 구분
- 6번, 7번, 8번: state equation 차원
- 27번, 28번, 30번: structural damping 적용 범위
- 41번, 42번, 44번: stable vs asymptotically stable
- 51번, 52번: proportional damping vs non-proportional damping
- 56번, 57번, 58번, 59번, 60번: gyroscopic / nonsymmetric EVP / bi-orthogonality

