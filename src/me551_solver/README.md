# ME551 Linear Vibration Midterm Solver

An interactive CLI solver for ME551 (Linear Vibration Engineering) midterm problems.
Runs fully offline — no internet connection required on exam day.

## How to Run

```bash
# From the src/ directory
cd src
python -m me551_solver
```

## How to Build a Standalone Executable

### Windows
```bat
cd src
build.bat
```

### Linux / macOS
```bash
cd src
bash build.sh
```

The built executable will be at `dist/me551_solver/me551_solver.exe` (Windows)
or `dist/me551_solver/me551_solver` (Linux/macOS).

You can also run PyInstaller directly:
```bash
pip install pyinstaller
pyinstaller me551_solver.spec --clean
```

## Available Problem Types (Main Menu)

| # | Problem Type |
|---|---|
| 1 | T/F 개념 판별 — keyword-based True/False verdict from lecture notes |
| 2 | FRF 분해 — transfer-function partial-fraction decomposition |
| 3 | Lagrange -> equilibrium -> linearization -> stability analysis |
| 4 | 2-DOF Modal Analysis (eigenvalues, mode shapes, free response) |
| 5 | Proportional Damping Response (Rayleigh damping, modal damping) |
| 6 | Gyroscopic / Nonconservative extension (v2) |
| 7 | Template: Rotating Hoop (bead on rotating hoop) |
| 8 | Template: Rotating Triangle (mass on rotating frame) |
| 9 | Template: 2-DOF Chain System (spring-mass chain) |

Press `h` at the menu for detailed help on each option.

## Dependencies

| Package | Version |
|---------|---------|
| sympy | >=1.12 |
| numpy | >=1.24 |
| scipy | >=1.10 |
| matplotlib | >=3.7 |
| pyyaml | >=6.0 |

Install with:
```bash
pip install -r requirements.txt
```

## Data Files

- `data/tf_rules.yaml` — True/False concept rules database (lecture-indexed)
- `data/lecture_index.yaml` — lecture topic index

These files are bundled automatically by PyInstaller via `me551_solver.spec`.
