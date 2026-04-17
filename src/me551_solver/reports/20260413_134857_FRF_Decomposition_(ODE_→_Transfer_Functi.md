# FRF Decomposition (ODE → Transfer Function → Partial Fractions)

## Given
- **ode_string**: x'''+1.2x''+4.2x'+4x=f(t)
- **coefficients**: [1.0, 1.2, 4.2, 4.0]
- **system_order**: 3
- **numerator**: [1]

## Derivation Steps
### Step 1: ODE string parsed

$$
Input:   'x'''+1.2x''+4.2x'+4x=f(t)'
Coefficients [a_n, …, a_0]: [1.0, 1.2, 4.2, 4.0]
$$

### Step 2: Transfer function G(s) 구성

$$
G(s) = 1 / (        2           
 3   6⋅s    21⋅s    
s  + ──── + ──── + 4
      5      5      )
$$

### Step 3: Characteristic equation roots (특성방정식의 근)

$$
Roots of         2           
 3   6⋅s    21⋅s    
s  + ──── + ──── + 4
      5      5       = 0:
  -1, -1/10 - sqrt(399)*I/10, -1/10 + sqrt(399)*I/10
$$

### Step 4: Partial fraction decomposition (부분분수 분해)

```
     5⋅(5⋅s - 4)           5     
- ────────────────── + ──────────
     ⎛   2         ⎞   24⋅(s + 1)
  24⋅⎝5⋅s  + s + 20⎠             
```

### Step 5: Denominator factorisation (분모 인수분해)

```
        ⎛   2         ⎞
(s + 1)⋅⎝5⋅s  + s + 20⎠
───────────────────────
           5           
```

### Step 6: 1st-order subsystem

$$
G₁(s) =     5     
──────────
24⋅(s + 1)
  Break frequency ω_b = 1
  Time constant τ = 1

  G₁(jω) = [5] / [(24) + j*(24*omega)]
  Re[G₁(jω)] = 5/(24*(omega**2 + 1))
  Im[G₁(jω)] = -5*omega/(24*omega**2 + 24)
$$

### Step 7: 2nd-order subsystem

$$
G₂(s) =   -5⋅(5⋅s - 4)    
──────────────────
   ⎛   2         ⎞
24⋅⎝5⋅s  + s + 20⎠
  Natural frequency ω_n = 2 (ω_n² = 4)
  Damping ratio ζ = 1/20

  G₂(jω) = [(20) + j*(-25*omega)] / [(480 - 120*omega**2) + j*(24*omega)]
  Re[G₂(jω)] = 25*(16 - 5*omega**2)/(24*(25*omega**4 - 199*omega**2 + 400))
  Im[G₂(jω)] = (125*omega**3 - 520*omega)/(600*omega**4 - 4776*omega**2 + 9600)
$$

### Step 8: DC gain |G(0)|

$$
|G(0)| = |1/4| = 1/4
$$

### Step 9: Characteristic points (특성점 요약)

$$
{'dc_gain': 0.25, 'first_order_subsystems': [{'omega_b': 1.0, 'mag_at_omega_b': 0.1473139127471974, 'description': 'ω_b = 1 rad/s, |G₁(jω_b)| = 0.1473'}], 'second_order_subsystems': [{'omega_n': 2.0, 'zeta': 0.05, 'mag_at_omega_n': 1.1219093348196882, 'description': 'ω_n = 2 rad/s, ζ = 0.05, |G₂(jω_n)| = 1.122'}], 'high_freq_rolloff_dB_per_decade': -60, 'high_freq_description': 'High-frequency asymptote: -60 dB/decade (relative order = 3)'}
$$

### Step 10: Summary

$$
System order: 3
Characteristic polynomial: 1/(s**3 + 6*s**2/5 + 21*s/5 + 4)
G₁(s): 1st-order,  ω_b = 1 rad/s
G₂(s): 2nd-order,  ω_n = 2 rad/s,  ζ = 0.05
$$

## Final Answer
- **transfer_function**: `1/(s**3 + 6*s**2/5 + 21*s/5 + 4)`
- **roots**: `['-1', '-1/10 - sqrt(399)*I/10', '-1/10 + sqrt(399)*I/10']`
- **partial_fractions**: `-5*(5*s - 4)/(24*(5*s**2 + s + 20)) + 5/(24*(s + 1))`
- **dc_gain**: `0.25`
- **characteristic_points**: `{'dc_gain': 0.25, 'first_order_subsystems': [{'omega_b': 1.0, 'mag_at_omega_b': 0.1473139127471974, 'description': 'ω_b = 1 rad/s, |G₁(jω_b)| = 0.1473'}], 'second_order_subsystems': [{'omega_n': 2.0, 'zeta': 0.05, 'mag_at_omega_n': 1.1219093348196882, 'description': 'ω_n = 2 rad/s, ζ = 0.05, |G₂(jω_n)| = 1.122'}], 'high_freq_rolloff_dB_per_decade': -60, 'high_freq_description': 'High-frequency asymptote: -60 dB/decade (relative order = 3)'}`
- **first_order**: `[{'omega_b': 1.0, 'time_constant': 1.0, 'subsystem_tf': '5/(24*(s + 1))', 'G_jw': '[5] / [(24) + j*(24*omega)]', 'G_jw_real': '5/(24*(omega**2 + 1))', 'G_jw_imag': '-5*omega/(24*omega**2 + 24)'}]`
- **second_order**: `[{'omega_n': 2.0, 'zeta': 0.05, 'subsystem_tf': '-5*(5*s - 4)/(24*(5*s**2 + s + 20))', 'G_jw': '[(20) + j*(-25*omega)] / [(480 - 120*omega**2) + j*(24*omega)]', 'G_jw_real': '25*(16 - 5*omega**2)/(24*(25*omega**4 - 199*omega**2 + 400))', 'G_jw_imag': '(125*omega**3 - 520*omega)/(600*omega**4 - 4776*omega**2 + 9600)'}]`
- **pole_real**: `-1.0`
- **bandwidth_G1**: `1.0`
- **omega_n**: `2.0`
- **zeta**: `0.05`
- **omega_1**: `2.0`
- **zeta_1**: `0.05`
- **G_dc**: `0.25`

## Sanity Check
Stability: All roots have negative real parts → STABLE
DC gain check: numerator_const/denominator_const = 1.0/4.0 = 0.25 ✓
Partial fraction recombination matches original G(s) ✓
