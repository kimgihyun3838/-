# FRF Decomposition (ODE → Transfer Function → Partial Fractions)

## Given
- **ode_string**: x'''+3x''+6x'+8x = f(t)
- **coefficients**: [1.0, 3.0, 6.0, 8.0]
- **system_order**: 3
- **numerator**: [1]

## Derivation Steps
### Step 1: ODE string parsed

$$
Input:   'x'''+3x''+6x'+8x = f(t)'
Coefficients [a_n, …, a_0]: [1.0, 3.0, 6.0, 8.0]
$$

### Step 2: Transfer function G(s) 구성

$$
G(s) = 1 / ( 3      2          
s  + 3⋅s  + 6⋅s + 8)
$$

### Step 3: Characteristic equation roots (특성방정식의 근)

$$
Roots of  3      2          
s  + 3⋅s  + 6⋅s + 8 = 0:
  -2, -1/2 - sqrt(15)*I/2, -1/2 + sqrt(15)*I/2
$$

### Step 4: Partial fraction decomposition (부분분수 분해)

```
      s - 1            1    
- ────────────── + ─────────
    ⎛ 2        ⎞   6⋅(s + 2)
  6⋅⎝s  + s + 4⎠            
```

### Step 5: Denominator factorisation (분모 인수분해)

```
        ⎛ 2        ⎞
(s + 2)⋅⎝s  + s + 4⎠
```

### Step 6: 1st-order subsystem

$$
G₁(s) =     1    
─────────
6⋅(s + 2)
  Break frequency ω_b = 2
  Time constant τ = 1/2

  G₁(jω) = [1] / [(12) + j*(6*omega)]
  Re[G₁(jω)] = 1/(3*(omega**2 + 4))
  Im[G₁(jω)] = -omega/(6*omega**2 + 24)
$$

### Step 7: 2nd-order subsystem

$$
G₂(s) =   -(s - 1)    
──────────────
  ⎛ 2        ⎞
6⋅⎝s  + s + 4⎠
  Natural frequency ω_n = 2 (ω_n² = 4)
  Damping ratio ζ = 1/4

  G₂(jω) = [(1) + j*(-omega)] / [(24 - 6*omega**2) + j*(6*omega)]
  Re[G₂(jω)] = (2 - omega**2)/(3*(omega**4 - 7*omega**2 + 16))
  Im[G₂(jω)] = omega*(omega**2 - 5)/(6*(omega**4 - 7*omega**2 + 16))
$$

### Step 8: DC gain |G(0)|

$$
|G(0)| = |1/8| = 1/8
$$

### Step 9: Characteristic points (특성점 요약)

$$
{'dc_gain': 0.125, 'first_order_subsystems': [{'omega_b': 2.0, 'mag_at_omega_b': 0.05892556509887896, 'description': 'ω_b = 2 rad/s, |G₁(jω_b)| = 0.05893'}], 'second_order_subsystems': [{'omega_n': 2.0, 'zeta': 0.25, 'mag_at_omega_n': 0.18633899812498247, 'description': 'ω_n = 2 rad/s, ζ = 0.25, |G₂(jω_n)| = 0.1863'}], 'high_freq_rolloff_dB_per_decade': -60, 'high_freq_description': 'High-frequency asymptote: -60 dB/decade (relative order = 3)'}
$$

### Step 10: Summary

$$
System order: 3
Characteristic polynomial: 1/(s**3 + 3*s**2 + 6*s + 8)
G₁(s): 1st-order,  ω_b = 2 rad/s
G₂(s): 2nd-order,  ω_n = 2 rad/s,  ζ = 0.25
$$

## Final Answer
- **transfer_function**: `1/(s**3 + 3*s**2 + 6*s + 8)`
- **roots**: `['-2', '-1/2 - sqrt(15)*I/2', '-1/2 + sqrt(15)*I/2']`
- **partial_fractions**: `-(s - 1)/(6*(s**2 + s + 4)) + 1/(6*(s + 2))`
- **dc_gain**: `0.125`
- **characteristic_points**: `{'dc_gain': 0.125, 'first_order_subsystems': [{'omega_b': 2.0, 'mag_at_omega_b': 0.05892556509887896, 'description': 'ω_b = 2 rad/s, |G₁(jω_b)| = 0.05893'}], 'second_order_subsystems': [{'omega_n': 2.0, 'zeta': 0.25, 'mag_at_omega_n': 0.18633899812498247, 'description': 'ω_n = 2 rad/s, ζ = 0.25, |G₂(jω_n)| = 0.1863'}], 'high_freq_rolloff_dB_per_decade': -60, 'high_freq_description': 'High-frequency asymptote: -60 dB/decade (relative order = 3)'}`
- **first_order**: `[{'omega_b': 2.0, 'time_constant': 0.5, 'subsystem_tf': '1/(6*(s + 2))', 'G_jw': '[1] / [(12) + j*(6*omega)]', 'G_jw_real': '1/(3*(omega**2 + 4))', 'G_jw_imag': '-omega/(6*omega**2 + 24)'}]`
- **second_order**: `[{'omega_n': 2.0, 'zeta': 0.25, 'subsystem_tf': '-(s - 1)/(6*(s**2 + s + 4))', 'G_jw': '[(1) + j*(-omega)] / [(24 - 6*omega**2) + j*(6*omega)]', 'G_jw_real': '(2 - omega**2)/(3*(omega**4 - 7*omega**2 + 16))', 'G_jw_imag': 'omega*(omega**2 - 5)/(6*(omega**4 - 7*omega**2 + 16))'}]`
- **pole_real**: `-2.0`
- **bandwidth_G1**: `2.0`
- **omega_n**: `2.0`
- **zeta**: `0.25`
- **omega_1**: `2.0`
- **zeta_1**: `0.25`
- **G_dc**: `0.125`

## Sanity Check
Stability: All roots have negative real parts → STABLE
DC gain check: numerator_const/denominator_const = 1.0/8.0 = 0.125 ✓
Partial fraction recombination matches original G(s) ✓
