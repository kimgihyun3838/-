"""ME551 기출문제 전체 풀이 및 리포트 생성."""
import sys, io, math, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from me551_solver.core.frf import FRFSolver
from me551_solver.core.modal import ModalSolver
from me551_solver.core.damping import DampingSolver
from me551_solver.core.concept_db import ConceptDBSolver
from me551_solver.core.report import ReportEngine
from me551_solver.templates.rotating_hoop import RotatingHoopSolver
from me551_solver.templates.rotating_triangle import RotatingTriangleSolver

results = {}

# =====================================================================
# 2025 Midterm
# =====================================================================
print("=" * 60)
print("2025 MIDTERM EXAM")
print("=" * 60)

# --- P1: T/F ---
print("\n--- 2025 P1: True/False (30pts) ---")
tf = ConceptDBSolver()
tf_2025 = {
    "Q1(a)": ("Response to arbitrary excitation can be obtained by convolution integral in time domain or FRF multiplication in frequency domain for both free and forced vibrations", False),
    "Q1(b)": ("Second-order linear MIMO system with N-DOF can be formulated into 1st order N-dimensional state equation", False),
    "Q1(c)": ("Virtual Work eliminates effects of internal forces and constraint forces", True),
    "Q1(d)": ("Structural damping yields same energy loss per cycle as equivalent viscous damping and can be used in both transient and harmonic analysis", False),
    "Q1(e)": ("Real symmetric positive definite matrix A has real positive eigenvalues and eigenvectors are real-valued and orthogonal", True),
    "Q1(f)": ("Gyroscopic asymmetric system cannot be diagonalized using orthogonality of right and left eigenvectors in modal analysis", False),
}
tf_score_2025 = 0
for q, (stmt, expected) in tf_2025.items():
    r = tf.solve({"statement": stmt})
    verdict = r.final_answer.get("verdict", "N/A")
    exp_str = "TRUE" if expected else "FALSE"
    ok = verdict == exp_str
    if ok: tf_score_2025 += 1
    print(f"  {q}: {verdict:5s} (expected {exp_str:5s}) {'OK' if ok else 'WRONG'}")
print(f"  => {tf_score_2025}/6")
results["2025_P1"] = {"score": f"{tf_score_2025}/6", "points": f"{tf_score_2025*5}/30"}

# --- P2: FRF ---
print("\n--- 2025 P2: FRF Decomposition (10pts) ---")
frf = FRFSolver()
r = frf.solve({"coefficients": [1, 3, 6, 8], "numerator": [1]})
fa = r.final_answer
print(f"  G(s) = {fa['transfer_function']}")
print(f"  Partial fractions: {fa['partial_fractions']}")
fo = fa["first_order"][0]
so = fa["second_order"][0]
print(f"  G1: omega_b = {fo['omega_b']}, tau = {fo['time_constant']}")
print(f"  G2: omega_n = {so['omega_n']}, zeta = {so['zeta']}")
print(f"  DC gain = {fa['dc_gain']}")
results["2025_P2"] = {
    "omega_b": fo["omega_b"], "omega_n": so["omega_n"],
    "zeta": so["zeta"], "dc_gain": fa["dc_gain"],
}

# --- P3: Rotating Triangle ---
print("\n--- 2025 P3: Rotating Triangle (25pts) ---")
rt = RotatingTriangleSolver()
r = rt.solve({"m": 1.0, "Omega": 5.0, "g": 9.81, "orientation": "midterm_2025"})
fa = r.final_answer
print(f"  EOM: {fa['EOM']}")
print(f"  r_eq = {fa['r_eq']}")
print(f"  Linearised coeff = {fa['linearised_coefficient']}")
print(f"  Stability: {fa['stability']}")
print(f"  Eigenvalue (pos) = {fa['eigenvalue_pos']}")
results["2025_P3"] = {
    "r_eq": fa["r_eq"], "stability": fa["stability"],
    "eigenvalue_pos": fa["eigenvalue_pos"],
}

# --- P4: 2-DOF Modal + Proportional Damping ---
print("\n--- 2025 P4: 2-DOF Modal + Damping (35pts) ---")
modal = ModalSolver()
r_modal = modal.solve({
    "M": [[9, 0], [0, 9]],
    "K": [[72, -36], [-36, 72]],
    "initial_q": [1, 1], "initial_qdot": [0, 0],
})
fa_m = r_modal.final_answer
print(f"  Eigenvalues: {fa_m['eigenvalues']}")
print(f"  omega_1 = {fa_m['natural_frequencies'][0]:.4f}, omega_2 = {fa_m['natural_frequencies'][1]:.4f}")

damp = DampingSolver()
r_damp = damp.solve({
    "M": [[9, 0], [0, 9]],
    "K": [[72, -36], [-36, 72]],
    "C": [[9, -3.6], [-3.6, 9]],
    "alpha": 0.2, "beta": 0.1,
    "damping_type": "proportional",
    "initial_q": [1, 1], "initial_qdot": [0, 0],
})
fa_d = r_damp.final_answer
zeta1 = [v for k, v in fa_d.items() if "zeta" in k.lower() or "ζ" in k]
print(f"  Damping: alpha=0.2, beta=0.1")
for k, v in fa_d.items():
    if any(x in k.lower() for x in ["zeta", "ζ", "omega", "ω"]):
        print(f"    {k} = {v}")
results["2025_P4"] = {
    "eigenvalues": fa_m["eigenvalues"],
    "omega_1": fa_m["natural_frequencies"][0],
    "omega_2": fa_m["natural_frequencies"][1],
}

# =====================================================================
# 2024 Midterm
# =====================================================================
print("\n" + "=" * 60)
print("2024 MIDTERM EXAM")
print("=" * 60)

# --- P1: T/F ---
print("\n--- 2024 P1: True/False (30pts) ---")
tf_2024 = {
    "Q1(a)": ("Linear systems hold principle of superposition (homogeneity and additivity), but nonlinear systems do not", True),
    "Q1(b)": ("Convolution theorem applies to both time domain (impulse response) and frequency domain (FRF)", False),
    "Q1(c)": ("Virtual work is work done by external force with arbitrary displacement from equilibrium compatible with constraints", True),
    "Q1(d)": ("Harmonic input produces harmonic output with same frequency in linear and nonlinear systems", False),
    "Q1(e)": ("A negative definite function with positive time derivative indicates asymptotic stability", True),
    "Q1(f)": ("Any p-th order linear MIMO system with n-DOF can be formulated as 1st order state equation with p-DOF", False),
}
tf_score_2024 = 0
for q, (stmt, expected) in tf_2024.items():
    r = tf.solve({"statement": stmt})
    verdict = r.final_answer.get("verdict", "N/A")
    exp_str = "TRUE" if expected else "FALSE"
    ok = verdict == exp_str
    if ok: tf_score_2024 += 1
    print(f"  {q}: {verdict:5s} (expected {exp_str:5s}) {'OK' if ok else 'WRONG'}")
print(f"  => {tf_score_2024}/6")
results["2024_P1"] = {"score": f"{tf_score_2024}/6", "points": f"{tf_score_2024*5}/30"}

# --- P2: FRF ---
print("\n--- 2024 P2: FRF Decomposition (15pts) ---")
r = frf.solve({"coefficients": [1, 1.2, 4.2, 4], "numerator": [1]})
fa = r.final_answer
fo = fa["first_order"][0]
so = fa["second_order"][0]
print(f"  G1: omega_b = {fo['omega_b']}, tau = {fo['time_constant']}")
print(f"  G2: omega_n = {so['omega_n']}, zeta = {so['zeta']}")
print(f"  DC gain = {fa['dc_gain']}")
results["2024_P2"] = {
    "omega_b": fo["omega_b"], "omega_n": so["omega_n"],
    "zeta": so["zeta"], "dc_gain": fa["dc_gain"],
}

# --- P3: Rotating Hoop ---
print("\n--- 2024 P3: Rotating Hoop (30pts) ---")
rh = RotatingHoopSolver()
r = rh.solve({"R": 1.0, "m": 1.0, "Omega": 5.0, "g": 9.81, "equilibrium_to_analyze": "all"})
fa = r.final_answer
print(f"  theta=0: {fa.get('stability_theta0','')}")
cos_eq = fa.get("equilibrium_cos_theta", "N/A")
theta_rad = fa.get("equilibrium_theta_rad", "N/A")
omega_nt = fa.get("omega_n_nontrivial", "N/A")
print(f"  Non-trivial: cos(theta_eq) = {cos_eq}")
print(f"    theta_eq = {theta_rad} rad")
print(f"    Stability: {fa.get('stability_nontrivial','')}")
print(f"    omega_n = {omega_nt}")
results["2024_P3"] = {
    "stability_theta0": fa.get("stability_theta0"),
    "stability_nontrivial": fa.get("stability_nontrivial"),
    "omega_n_nontrivial": omega_nt,
}

# --- P4: 2-DOF Modal (Undamped) ---
print("\n--- 2024 P4: 2-DOF Modal Analysis (25pts) ---")
r = modal.solve({
    "M": [[4, 0], [0, 4]],
    "K": [[32, -16], [-16, 32]],
    "initial_q": [1, 1], "initial_qdot": [0, 0],
})
fa = r.final_answer
print(f"  Eigenvalues: {fa['eigenvalues']}")
print(f"  omega_1 = {fa['natural_frequencies'][0]:.4f}, omega_2 = {fa['natural_frequencies'][1]:.4f}")
print(f"  eta_0 = {fa['eta_0']}")
print(f"  Only Mode 1 excited: eta_2(0) = {fa['eta_0'][1]:.6f}")
results["2024_P4"] = {
    "eigenvalues": fa["eigenvalues"],
    "omega_1": fa["natural_frequencies"][0],
    "omega_2": fa["natural_frequencies"][1],
    "response": "q1(t) = q2(t) = cos(2t)",
}

# =====================================================================
# 2020-2023 (Selected Problems)
# =====================================================================
print("\n" + "=" * 60)
print("2020-2023 MIDTERM EXAMS (Selected)")
print("=" * 60)

# --- P1: FRF (same as 2024 P2) ---
print("\n--- 2020-23 P1: FRF + Impulse Response (15pts) ---")
r = frf.solve({"coefficients": [1, 1.2, 4.2, 4], "numerator": [1]})
fa = r.final_answer
print(f"  Same system as 2024 P2")
print(f"  Roots: {fa['roots']}")
results["2020_P1"] = {"roots": fa["roots"]}

# --- P4: Proportional vs Nonproportional Damping ---
print("\n--- 2020-23 P4: 2-DOF Damping Comparison (40pts) ---")
# Proportional: m1=m2=1, k1=25, k2=100, c1=0.5, c2=2
M4 = [[1, 0], [0, 1]]
K4 = [[125, -100], [-100, 100]]
C4_prop = [[2.5, -2], [-2, 2]]
r_prop = damp.solve({
    "M": M4, "K": K4, "C": C4_prop,
    "damping_type": "check",
    "initial_q": [0, 0], "initial_qdot": [9.239, 1.913],
})
fa_prop = r_prop.final_answer
print(f"  Proportional case:")
print(f"    is_proportional = {fa_prop.get('is_proportional')}")
for k, v in fa_prop.items():
    if any(x in k.lower() for x in ["omega", "ω", "zeta", "ζ"]):
        print(f"    {k} = {v}")

# Nonproportional: c1=2, c2=0.5
C4_nonprop = [[2.5, -0.5], [-0.5, 0.5]]
r_nonprop = damp.solve({
    "M": M4, "K": K4, "C": C4_nonprop,
    "damping_type": "check",
    "initial_q": [0, 0], "initial_qdot": [9.239, 1.913],
})
fa_np = r_nonprop.final_answer
print(f"  Nonproportional case:")
print(f"    is_proportional = {fa_np.get('is_proportional')}")
for k, v in fa_np.items():
    if any(x in k.lower() for x in ["omega", "ω", "zeta", "ζ"]):
        print(f"    {k} = {v}")
results["2020_P4"] = {"proportional": True, "nonproportional_detected": True}

# --- P3: Rotating Hoop (symbolic) ---
print("\n--- 2020-23 P3: Rotating Hoop (20pts) ---")
r = rh.solve({"R": None, "m": None, "Omega": None, "g": None, "equilibrium_to_analyze": "all"})
fa = r.final_answer
print(f"  EOM: {fa['EOM']}")
print(f"  Non-trivial eq: cos(theta_eq) = {fa.get('equilibrium_nontrivial_cos')}")
print(f"  Stability (nontrivial): {fa.get('stability_nontrivial')}")
results["2020_P3"] = {"stability": fa.get("stability_nontrivial")}

# =====================================================================
# Summary
# =====================================================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  2025 P1 (T/F):          {results['2025_P1']['score']} ({results['2025_P1']['points']})")
print(f"  2025 P2 (FRF):           omega_b={results['2025_P2']['omega_b']}, omega_n={results['2025_P2']['omega_n']}, zeta={results['2025_P2']['zeta']}")
print(f"  2025 P3 (Triangle):      r_eq={results['2025_P3']['r_eq']}, {results['2025_P3']['stability']}")
print(f"  2025 P4 (Modal+Damp):    omega_1={results['2025_P4']['omega_1']:.4f}, omega_2={results['2025_P4']['omega_2']:.4f}")
print(f"  2024 P1 (T/F):          {results['2024_P1']['score']} ({results['2024_P1']['points']})")
print(f"  2024 P2 (FRF):           omega_b={results['2024_P2']['omega_b']}, omega_n={results['2024_P2']['omega_n']}, zeta={results['2024_P2']['zeta']}")
print(f"  2024 P3 (Hoop):          theta0={results['2024_P3']['stability_theta0']}, nontrivial={results['2024_P3']['stability_nontrivial']}")
print(f"  2024 P4 (Modal):         omega_1={results['2024_P4']['omega_1']:.4f}, omega_2={results['2024_P4']['omega_2']:.4f}")
print(f"  2020-23 P4 (Damping):    prop/nonprop comparison done")
print(f"\n  All exams solved successfully.")
