# Rotating Hoop – Lagrange / Equilibrium / Stability

## Given
- **R (hoop radius)**: symbolic R
- **m (bead mass)**: symbolic m
- **Omega (angular vel.)**: symbolic Omega
- **g (gravity)**: 9.81
- **Equilibria analysed**: all

## Derivation Steps
### Step 1: Kinetic energy T

$$
T = (1/2)*m*(R²·θ̇² + R²·Ω²·sin²θ)
  = m*(Omega**2*R**2*sin(theta)**2 + R**2*thetadot**2)/2
$$

### Step 2: Potential energy V

$$
V = -m·g·R·cos(θ)
  = -981*R*m*cos(theta)/100
$$

### Step 3: Lagrange's equation (EOM) — divide by m·R²

$$
θ'' + (g/R)·sin(θ) − Ω²·sin(θ)·cos(θ) = 0
θ'' + sin(θ)·[(g/R) − Ω²·cos(θ)] = 0
Symbolic form:
  theta'' + -Omega**2*sin(theta)*cos(theta) + 981*sin(theta)/(100*R) = 0
$$

### Step 4: Equilibrium conditions: set θ'' = 0 → sin(θ)·[(g/R) − Ω²cos(θ)] = 0

$$
Case A: sin(θ_eq) = 0  → θ_eq = 0  (bottom)
Case B: cos(θ_eq) = g/(R·Ω²)  — non-trivial (exists if |g/(R·Ω²)| ≤ 1)
$$

### Step 5: Non-trivial equilibrium (symbolic)

$$
cos(θ_eq) = g/(R·Ω²) = 981/(100*Omega**2*R)
Exists when Ω² > g/R  (i.e., Ω > √(g/R))
$$

### Step 6: Linearisation about θ_eq = 0

$$
Let θ = 0 + ε  (ε small)
sin(ε) ≈ ε,  cos(ε) ≈ 1
EOM → ε'' + [(g/R) − Ω²]·ε = 0
$$

### Step 7: Linearised coefficient at θ=0

$$
κ₀ = g/R − Ω²  =  -Omega**2 + 981/(100*R)
  κ₀ > 0  →  STABLE   (Ω < √(g/R))
  κ₀ < 0  →  UNSTABLE (Ω > √(g/R))
  κ₀ = 0  →  MARGINALLY STABLE (Ω = √(g/R))
$$

### Step 8: Critical speed for θ=0 stability

$$
Ω_crit = √(g/R)  →  STABLE for Ω < Ω_crit
$$

### Step 9: Linearisation about θ_eq = arccos(g/(R·Ω²))

$$
Let θ = θ_eq + ε
Expand EOM to first order in ε:
  θ'' → ε''
  sin(θ_eq+ε)·[(g/R) − Ω²cos(θ_eq+ε)]
  ≈ [sin(θ_eq) + ε·cos(θ_eq)]·[(g/R) − Ω²cos(θ_eq) + ε·Ω²sin(θ_eq)]
At equilibrium: (g/R) − Ω²·cos(θ_eq) = 0
Collecting ε terms:
  ε'' + Ω²·sin²(θ_eq)·ε = 0
$$

### Step 10: Natural frequency at non-trivial equilibrium

$$
ω_n = Ω·|sin(θ_eq)|
    = Ω·√(1 − [g/(R·Ω²)]²)
    = sqrt(10000*Omega**4*R**2 - 962361)/(100*Omega*R)
Restoring coefficient: Ω²·sin²(θ_eq) = Omega**2 - 962361/(10000*Omega**2*R**2)
$$

## Final Answer
- **EOM**: `theta'' + (g/R)*sin(theta) - Omega^2*sin(theta)*cos(theta) = 0`
- **equilibrium_theta0**: `theta = 0 (always)`
- **equilibrium_nontrivial_cos**: `981/(100*Omega**2*R)`
- **stability_theta0**: `STABLE if Ω < √(g/R); UNSTABLE if Ω > √(g/R)`
- **linearised_coefficient_theta0**: `-Omega**2 + 981/(100*R)`
- **stability_nontrivial**: `STABLE (for Ω > √(g/R))`
- **omega_n_nontrivial**: `Omega*sqrt(1-(g/(R*Omega^2))^2)`
- **linearised_coefficient_nontrivial**: `Omega**2 - 962361/(10000*Omega**2*R**2)`

## Sanity Check
EOM verified: θ=0 satisfies sin(0)=0 → EOM vanishes ✓
