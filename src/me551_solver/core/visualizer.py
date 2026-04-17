"""ME551 시험용 시각화 모듈.

SolverResult 객체를 받아 교재 스타일의 공학 플롯을 생성한다.
matplotlib만 사용하며, 오프라인 환경에서 동작한다.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# matplotlib import with graceful degradation
# ---------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("TkAgg")          # interactive window; fallback handled below
except Exception:
    pass

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
    _MPL_AVAILABLE = True
except ImportError:
    _MPL_AVAILABLE = False

from .base import SolverResult

# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------
_STABLE_COLOR   = "#2563EB"   # blue
_UNSTABLE_COLOR = "#DC2626"   # red
_NEUTRAL_COLOR  = "#059669"   # green
_ACCENT_COLOR   = "#7C3AED"   # purple
_FIG_SINGLE     = (10, 6)
_FIG_SIDE       = (12, 5)


def _require_mpl() -> None:
    if not _MPL_AVAILABLE:
        raise ImportError(
            "matplotlib이 설치되어 있지 않습니다.  "
            "pip install matplotlib 으로 설치하세요."
        )


def _apply_style() -> None:
    """일관된 플롯 스타일을 적용한다."""
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        try:
            plt.style.use("seaborn-whitegrid")
        except OSError:
            pass  # 기본 스타일 사용

    plt.rcParams.update({
        "figure.dpi":        100,
        "axes.titlesize":    13,
        "axes.labelsize":    11,
        "legend.fontsize":   10,
        "xtick.labelsize":   9,
        "ytick.labelsize":   9,
        "lines.linewidth":   2.0,
        "grid.alpha":        0.4,
    })


def _safe_float(val: Any) -> float | None:
    """Safely convert any value to float, return None on failure."""
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def _parse_complex(s: str) -> complex | None:
    """Parse a sympy-style complex string to a Python complex number."""
    if s is None:
        return None
    s = str(s).strip()
    # Remove sympy sqrt, I, etc.
    try:
        # Replace sympy I with j, handle simple cases
        s2 = (s.replace("I", "j")
               .replace("sqrt", "cmath.sqrt")
               .replace("**", "**"))
        import cmath  # noqa: F401
        return complex(eval(s2, {"cmath": cmath, "__builtins__": {}}))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class Visualizer:
    """ME551 시험용 시각화 도구.

    Parameters
    ----------
    save_dir : str or Path, optional
        플롯 저장 디렉터리. 기본값은 ``reports/figures/``.
    """

    def __init__(self, save_dir: str | Path | None = None) -> None:
        _require_mpl()
        _apply_style()

        if save_dir is not None:
            self._save_dir = Path(save_dir)
        else:
            # Resolve relative to this file: …/me551_solver/reports/figures/
            self._save_dir = (
                Path(__file__).resolve().parent.parent / "reports" / "figures"
            )

    # ------------------------------------------------------------------
    # Internal save helper
    # ------------------------------------------------------------------

    def _save_fig(self, fig: "plt.Figure", name: str) -> None:
        """Save figure to save_dir with the given base name."""
        self._save_dir.mkdir(parents=True, exist_ok=True)
        path = self._save_dir / f"{name}.png"
        fig.savefig(path, bbox_inches="tight", dpi=150)
        print(f"  [저장됨] {path}")

    # ------------------------------------------------------------------
    # 1. Bode Plot
    # ------------------------------------------------------------------

    def plot_bode(
        self,
        frf_result: SolverResult,
        save: bool = False,
    ) -> "plt.Figure":
        """Bode 선도 (Magnitude dB + Phase deg) 를 그린다.

        Parameters
        ----------
        frf_result : SolverResult
            FRFSolver 또는 FRFDecomposeSolver 의 출력.
        save : bool
            True 이면 save_dir 에 PNG 저장.

        Returns
        -------
        matplotlib.figure.Figure
        """
        _apply_style()
        fa = frf_result.final_answer
        given = frf_result.given

        # ---- Frequency data ----
        omega = fa.get("omega")
        mag_dB = fa.get("magnitude_dB")
        phase_deg = fa.get("phase_deg")

        # If no pre-computed frequency sweep, build one from coefficients
        if omega is None or mag_dB is None:
            den = given.get("denominator_coefficients") or given.get("coefficients")
            num = given.get("numerator_coefficients") or given.get("numerator", [1])
            if den is None:
                raise ValueError(
                    "BodePlot: FRF result에 omega/magnitude_dB 데이터가 없고 "
                    "계수 정보도 없습니다. FRFSolver에 omega_range를 지정하거나 "
                    "FRFDecomposeSolver를 사용하세요."
                )
            den_f = [float(c) for c in den]
            num_f = [float(c) for c in num]

            # Guess a sensible range from characteristic points
            cp = fa.get("characteristic_points", {})
            w_pts = []
            for s1 in cp.get("first_order_subsystems", []):
                if "omega_b" in s1:
                    w_pts.append(float(s1["omega_b"]))
            for s2 in cp.get("second_order_subsystems", []):
                if "omega_n" in s2:
                    w_pts.append(float(s2["omega_n"]))
            if w_pts:
                w_min = min(w_pts) / 100.0
                w_max = max(w_pts) * 100.0
            else:
                w_min, w_max = 0.01, 1000.0

            omega = np.logspace(np.log10(w_min), np.log10(w_max), 600)
            jw = 1j * omega
            num_v = np.polyval([complex(c) for c in num_f], jw)
            den_v = np.polyval([complex(c) for c in den_f], jw)
            G_jw = num_v / den_v
            mag_dB   = 20 * np.log10(np.abs(G_jw))
            phase_deg = np.degrees(np.angle(G_jw))
        else:
            omega     = np.asarray(omega)
            mag_dB    = np.asarray(mag_dB)
            phase_deg = np.asarray(phase_deg)

        # ---- Build figure ----
        fig, (ax_mag, ax_ph) = plt.subplots(
            2, 1, figsize=_FIG_SINGLE, sharex=True,
            gridspec_kw={"hspace": 0.05},
        )

        tf_str = fa.get("transfer_function", "G(s)")
        fig.suptitle(f"Bode Plot — G(s) = {tf_str}", fontsize=12, y=1.01)

        # ---- Magnitude (linear scale) ----
        mag_linear = 10 ** (mag_dB / 20.0)
        ax_mag.loglog(omega, mag_linear, color=_STABLE_COLOR, label="G(j\u03c9)", zorder=3)

        # DC gain line
        dc_gain = _safe_float(fa.get("dc_gain"))
        if dc_gain is not None and dc_gain > 0:
            ax_mag.axhline(dc_gain, color="gray", linestyle=":", linewidth=1.2,
                           label=f"DC gain = {dc_gain:.3g}")

        # Characteristic points
        cp = fa.get("characteristic_points", {})
        colors_cp = [_ACCENT_COLOR, _UNSTABLE_COLOR, _NEUTRAL_COLOR]

        for i, s1 in enumerate(cp.get("first_order_subsystems", [])):
            wb = _safe_float(s1.get("omega_b"))
            if wb and wb > 0:
                clr = colors_cp[i % len(colors_cp)]
                ax_mag.axvline(wb, color=clr, linestyle="--", linewidth=1.2, alpha=0.8,
                               label=f"\u03c9_b = {wb:.3g} rad/s")

        for i, s2 in enumerate(cp.get("second_order_subsystems", [])):
            wn = _safe_float(s2.get("omega_n"))
            if wn and wn > 0:
                clr = colors_cp[(i + 1) % len(colors_cp)]
                ax_mag.axvline(wn, color=clr, linestyle="--", linewidth=1.2, alpha=0.8,
                               label=f"\u03c9_n = {wn:.3g} rad/s")

        # Subsystem contributions (dashed)
        self._plot_subsystem_contributions(ax_mag, fa, omega)

        ax_mag.set_ylabel("Magnitude |G(j\u03c9)|")
        ax_mag.legend(loc="lower left", fontsize=8)
        ax_mag.grid(True, which="both", alpha=0.3)
        ax_mag.set_title("Magnitude Response")

        # ---- Phase ----
        ax_ph.semilogx(omega, phase_deg, color=_ACCENT_COLOR, label="Phase")
        ax_ph.axhline(-90, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
        ax_ph.axhline(-180, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
        ax_ph.set_xlabel("\u03c9 (rad/s)")
        ax_ph.set_ylabel("Phase (deg)")
        ax_ph.legend(loc="lower left", fontsize=8)
        ax_ph.grid(True, which="both", alpha=0.3)
        ax_ph.set_title("Phase Response")

        fig.tight_layout()

        if save:
            self._save_fig(fig, "bode_plot")

        return fig

    def _plot_subsystem_contributions(
        self,
        ax: "plt.Axes",
        fa: dict,
        omega: np.ndarray,
    ) -> None:
        """1차/2차 서브시스템의 개별 G(j\u03c9) magnitude를 파선으로 그린다."""
        styles = ["--", "-."]
        colors = ["#F59E0B", "#10B981"]

        for i, s1 in enumerate(fa.get("first_order", [])):
            wb = _safe_float(s1.get("omega_b"))
            if wb is None or wb == 0:
                continue
            tf_str = s1.get("subsystem_tf", "")
            g1 = self._eval_tf_string(tf_str, omega)
            if g1 is not None:
                mag_lin = np.abs(g1) + 1e-300
                ax.loglog(omega, mag_lin,
                          linestyle=styles[i % len(styles)],
                          color=colors[i % len(colors)],
                          linewidth=1.2, alpha=0.75,
                          label=f"G{i+1}(j\u03c9)")

        for i, s2 in enumerate(fa.get("second_order", [])):
            tf_str = s2.get("subsystem_tf", "")
            g2 = self._eval_tf_string(tf_str, omega)
            if g2 is not None:
                mag_lin = np.abs(g2) + 1e-300
                ax.loglog(omega, mag_lin,
                          linestyle=styles[(i + 1) % len(styles)],
                          color=colors[(i + 1) % len(colors)],
                          linewidth=1.2, alpha=0.75,
                          label=f"G{i+2}(j\u03c9)")

    @staticmethod
    def _eval_tf_string(tf_str: str, omega: np.ndarray) -> np.ndarray | None:
        """Attempt to evaluate a sympy-format TF string G(s) at s=jw."""
        if not tf_str:
            return None
        try:
            import sympy as sp
            s_sym = sp.Symbol("s")
            expr = sp.sympify(tf_str)
            f = sp.lambdify(s_sym, expr, modules="numpy")
            return f(1j * omega)
        except Exception:
            return None

    # ------------------------------------------------------------------
    # 2. Pole-Zero Map
    # ------------------------------------------------------------------

    def plot_pole_zero(
        self,
        result: SolverResult,
        save: bool = False,
    ) -> "plt.Figure":
        """극점-영점 지도(Pole-Zero Map)를 그린다.

        Parameters
        ----------
        result : SolverResult
            final_answer에 'roots', 'eigenvalues', 또는 'eigenvalues_numeric'
            키가 있는 모든 SolverResult.
        save : bool
            True 이면 save_dir 에 PNG 저장.
        """
        _apply_style()
        fa = result.final_answer

        # Collect poles (roots of denominator)
        raw_poles = (
            fa.get("roots")
            or fa.get("eigenvalues_numeric")
            or fa.get("eigenvalues")
        )
        if raw_poles is None:
            raise ValueError(
                "plot_pole_zero: result.final_answer에 'roots', 'eigenvalues', "
                "또는 'eigenvalues_numeric' 키가 없습니다."
            )

        poles: list[complex] = []
        for p in raw_poles:
            if isinstance(p, complex):
                poles.append(p)
            else:
                c = _parse_complex(str(p))
                if c is not None:
                    poles.append(c)

        if not poles:
            raise ValueError("plot_pole_zero: 유효한 극점 데이터가 없습니다.")

        # Collect zeros if available
        zeros_raw = fa.get("zeros", [])
        zeros: list[complex] = []
        for z in zeros_raw:
            if isinstance(z, complex):
                zeros.append(z)
            else:
                c = _parse_complex(str(z))
                if c is not None:
                    zeros.append(c)

        # ---- Figure ----
        fig, ax = plt.subplots(figsize=_FIG_SINGLE)
        fig.suptitle("Pole-Zero Map", fontsize=13)

        # ---- Background: stable (left) vs unstable (right) ----
        xlim_margin = max(
            max(abs(p.real) for p in poles) * 2.5,
            max((abs(p.imag) for p in poles), default=1.0) * 0.5,
            1.0,
        )
        ylim_margin = max(
            max(abs(p.imag) for p in poles) * 2.5,
            1.0,
        )

        ax.axvspan(-xlim_margin, 0, alpha=0.06, color=_NEUTRAL_COLOR, zorder=0)
        ax.axvspan(0, xlim_margin, alpha=0.06, color=_UNSTABLE_COLOR, zorder=0)

        # ---- Axes ----
        ax.axhline(0, color="black", linewidth=0.8, zorder=1)
        ax.axvline(0, color="black", linewidth=1.2, zorder=1, label="Im axis")

        # ---- omega_n circles ----
        omega_n_vals = sorted(set(round(abs(p), 4) for p in poles if abs(p) > 1e-6))
        for wn in omega_n_vals:
            theta = np.linspace(0, 2 * np.pi, 200)
            ax.plot(wn * np.cos(theta), wn * np.sin(theta),
                    color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
            ax.annotate(f"\u03c9_n={wn:.3g}",
                        xy=(wn * np.cos(np.pi * 0.15), wn * np.sin(np.pi * 0.15)),
                        fontsize=7, color="gray", ha="left")

        # ---- zeta lines (constant damping ratio) ----
        zeta_vals = []
        for p in poles:
            wn = abs(p)
            if wn > 1e-6 and p.real < -1e-6:
                z = -p.real / wn
                if 0 < z < 1:
                    zeta_vals.append(round(z, 3))
        for zv in sorted(set(zeta_vals)):
            angle = np.arccos(zv)  # angle from negative real axis
            r_line = xlim_margin
            ax.plot(
                [-r_line * np.cos(angle), 0],
                [r_line * np.sin(angle), 0],
                color="gray", linestyle=":", linewidth=0.8, alpha=0.5,
            )
            ax.plot(
                [-r_line * np.cos(angle), 0],
                [-r_line * np.sin(angle), 0],
                color="gray", linestyle=":", linewidth=0.8, alpha=0.5,
            )
            ax.annotate(f"\u03b6={zv:.2f}",
                        xy=(-r_line * 0.7 * np.cos(angle),
                            r_line * 0.7 * np.sin(angle)),
                        fontsize=7, color="gray")

        # ---- Plot poles ----
        for p in poles:
            color = _STABLE_COLOR if p.real <= 1e-8 else _UNSTABLE_COLOR
            ax.plot(p.real, p.imag, "x", color=color,
                    markersize=12, markeredgewidth=2.5, zorder=5)
            ax.annotate(
                f" ({p.real:.3g}{'+' if p.imag >= 0 else ''}{p.imag:.3g}j)",
                xy=(p.real, p.imag),
                fontsize=7.5, color=color, zorder=6,
            )

        # ---- Plot zeros ----
        for z in zeros:
            ax.plot(z.real, z.imag, "o", color=_ACCENT_COLOR,
                    markersize=9, markerfacecolor="none",
                    markeredgewidth=2, zorder=5)
            ax.annotate(
                f" ({z.real:.3g}{'+' if z.imag >= 0 else ''}{z.imag:.3g}j)",
                xy=(z.real, z.imag),
                fontsize=7.5, color=_ACCENT_COLOR, zorder=6,
            )

        # ---- Legend proxies ----
        pole_stable_patch = plt.Line2D([], [], marker="x", color=_STABLE_COLOR,
                                       linestyle="None", markersize=10,
                                       markeredgewidth=2.5, label="Stable pole")
        pole_unstable_patch = plt.Line2D([], [], marker="x", color=_UNSTABLE_COLOR,
                                         linestyle="None", markersize=10,
                                         markeredgewidth=2.5, label="Unstable pole")
        zero_patch = plt.Line2D([], [], marker="o", color=_ACCENT_COLOR,
                                linestyle="None", markersize=9,
                                markerfacecolor="none", label="Zero")
        handles = [pole_stable_patch, pole_unstable_patch]
        if zeros:
            handles.append(zero_patch)
        ax.legend(handles=handles, loc="upper right", fontsize=9)

        ax.set_xlabel("Real part  \u03c3")
        ax.set_ylabel("Imaginary part  j\u03c9")
        ax.set_xlim(-xlim_margin, xlim_margin)
        ax.set_ylim(-ylim_margin, ylim_margin)
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal", adjustable="datalim")

        fig.tight_layout()

        if save:
            self._save_fig(fig, "pole_zero_map")

        return fig

    # ------------------------------------------------------------------
    # 3. Mode Shapes
    # ------------------------------------------------------------------

    def plot_mode_shapes(
        self,
        modal_result: SolverResult,
        save: bool = False,
    ) -> "plt.Figure":
        """2-DOF 이상 시스템의 모드형상을 막대 차트로 그린다.

        Parameters
        ----------
        modal_result : SolverResult
            ModalSolver 의 출력.  final_answer에 mass_normalized_eigenvectors,
            natural_frequencies 키가 필요하다.
        save : bool
            True 이면 save_dir 에 PNG 저장.
        """
        _apply_style()
        fa = modal_result.final_answer

        U_list = fa.get("mass_normalized_eigenvectors") or fa.get("modal_matrix")
        omega_list = fa.get("natural_frequencies")

        if U_list is None or omega_list is None:
            raise ValueError(
                "plot_mode_shapes: final_answer에 "
                "'mass_normalized_eigenvectors' 와 'natural_frequencies' 가 필요합니다."
            )

        U = np.array(U_list)
        omega = np.array(omega_list)
        n_dof, n_modes = U.shape

        # Normalise each mode to max |displacement| = 1
        U_norm = U.copy()
        for r in range(n_modes):
            peak = np.max(np.abs(U_norm[:, r]))
            if peak > 1e-12:
                U_norm[:, r] /= peak

        # ---- Figure: side-by-side ----
        fig, axes = plt.subplots(1, n_modes, figsize=_FIG_SIDE, sharey=True)
        if n_modes == 1:
            axes = [axes]

        fig.suptitle("Mode Shapes  (normalised to max |u| = 1)", fontsize=13)

        x_positions = np.arange(1, n_dof + 1)   # mass index
        colors_mode = [_STABLE_COLOR, _ACCENT_COLOR, _NEUTRAL_COLOR, _UNSTABLE_COLOR]

        for r, ax in enumerate(axes):
            shape = U_norm[:, r]
            color = colors_mode[r % len(colors_mode)]

            # Bar chart
            ax.barh(x_positions, shape, color=color, alpha=0.75,
                    edgecolor="black", linewidth=0.8)

            # Zero line
            ax.axvline(0, color="black", linewidth=1.0)

            # Connecting line (as in textbooks)
            ax.plot(shape, x_positions, "o-", color=color,
                    markersize=6, linewidth=1.5, zorder=3)

            # Label each bar
            for i, val in enumerate(shape):
                ax.text(val + 0.03 * np.sign(val) if abs(val) > 0.05 else 0.05,
                        x_positions[i], f"{val:.3f}",
                        va="center", ha="left" if val >= 0 else "right",
                        fontsize=8)

            freq_hz = omega[r] / (2 * np.pi)
            ax.set_title(
                f"Mode {r+1}\n"
                f"\u03c9\u2081 = {omega[r]:.4g} rad/s\n"
                f"({freq_hz:.3g} Hz)",
                fontsize=9,
            )
            ax.set_xlabel("Displacement (norm.)")
            ax.set_xlim(-1.5, 1.5)
            ax.set_yticks(x_positions)
            ax.set_yticklabels([f"DOF {i}" for i in x_positions])
            ax.grid(True, axis="x", alpha=0.3)

        fig.tight_layout()

        if save:
            self._save_fig(fig, "mode_shapes")

        return fig

    # ------------------------------------------------------------------
    # 4. Time Response
    # ------------------------------------------------------------------

    def plot_time_response(
        self,
        result: SolverResult,
        t_span: tuple[float, float] = (0.0, 5.0),
        n_points: int = 500,
        save: bool = False,
    ) -> "plt.Figure":
        """자유 진동 시간 응답 q(t) 를 그린다.

        ModalSolver 또는 DampingSolver 의 SolverResult를 받아
        각 DOF 의 시간 응답을 수치 계산하여 플롯한다.

        Parameters
        ----------
        result : SolverResult
        t_span : (t_start, t_end)
        n_points : int
        save : bool
        """
        _apply_style()
        fa = result.final_answer
        given = result.given

        t = np.linspace(float(t_span[0]), float(t_span[1]), n_points)
        is_damped = "Damped" in result.problem_type or "damped" in result.problem_type.lower()

        # ----- Determine DOF count -----
        n_dof = int(given.get("n_dof", 0))
        if n_dof == 0:
            # Infer from natural_frequencies
            nf = fa.get("natural_frequencies") or fa.get("undamped_natural_frequencies")
            if nf:
                n_dof = len(nf)
        if n_dof == 0:
            raise ValueError(
                "plot_time_response: DOF 수를 결정할 수 없습니다. "
                "ModalSolver 또는 DampingSolver 결과를 입력하세요."
            )

        # ----- Reconstruct response numerically -----
        U_list = fa.get("modal_matrix") or fa.get("mass_normalized_eigenvectors")
        if U_list is None:
            raise ValueError(
                "plot_time_response: final_answer에 'modal_matrix' 가 없습니다."
            )
        U = np.array(U_list)

        eta0 = np.array(fa.get("eta_0", np.zeros(n_dof)))
        etadot0 = np.array(fa.get("etadot_0", np.zeros(n_dof)))

        # Retrieve undamped omega
        omega_raw = fa.get("natural_frequencies") or fa.get("undamped_natural_frequencies")
        omega = np.array(omega_raw) if omega_raw else np.zeros(n_dof)

        # Retrieve damping ratios if available
        zeta = np.zeros(n_dof)
        omega_d = omega.copy()
        for r in range(n_dof):
            z = fa.get(f"zeta_{r+1}")
            if z is not None:
                zeta[r] = float(z)
            wd = fa.get(f"omega_d_{r+1}")
            if wd is not None:
                omega_d[r] = float(wd)
            # If omega_d not stored, compute it
            if fa.get(f"omega_d_{r+1}") is None:
                if zeta[r] < 1.0 and omega[r] > 1e-14:
                    omega_d[r] = omega[r] * np.sqrt(max(0, 1 - zeta[r] ** 2))

        # ---- Compute modal responses ----
        eta = np.zeros((n_dof, n_points))
        for r in range(n_dof):
            sigma_r = zeta[r] * omega[r]
            wd_r = omega_d[r]
            if omega[r] < 1e-14:
                # Rigid-body
                eta[r] = eta0[r] + etadot0[r] * t
            elif wd_r > 1e-14:
                # Underdamped
                A_r = eta0[r]
                B_r = (etadot0[r] + sigma_r * eta0[r]) / wd_r
                eta[r] = (
                    np.exp(-sigma_r * t)
                    * (A_r * np.cos(wd_r * t) + B_r * np.sin(wd_r * t))
                )
            else:
                # Undamped (no damping)
                A_r = eta0[r]
                B_r = etadot0[r] / omega[r] if omega[r] > 1e-14 else 0.0
                eta[r] = A_r * np.cos(omega[r] * t) + B_r * np.sin(omega[r] * t)

        # Physical response: q = U @ eta
        q = U @ eta   # shape (n_dof, n_points)

        # ---- Figure ----
        fig, ax = plt.subplots(figsize=_FIG_SINGLE)

        cmap_colors = [_STABLE_COLOR, _ACCENT_COLOR, _NEUTRAL_COLOR, _UNSTABLE_COLOR]
        for i in range(n_dof):
            color = cmap_colors[i % len(cmap_colors)]
            ax.plot(t, q[i], color=color, linewidth=2, label=f"q\u2081{i+1}(t)" if i == 0 else f"q_{i+1}(t)")

            # Damped envelope
            if is_damped and np.any(zeta > 1e-6):
                # Use modal contribution envelope per DOF (approximate)
                sigma_r = zeta[0] * omega[0] if omega[0] > 1e-14 else 0.0
                amp_approx = np.max(np.abs(q[i][:10]))
                if amp_approx > 1e-12 and sigma_r > 0:
                    envelope = amp_approx * np.exp(-sigma_r * t)
                    ax.plot(t, envelope, color=color, linestyle="--",
                            linewidth=1.0, alpha=0.5, label=f"envelope (DOF {i+1})")
                    ax.plot(t, -envelope, color=color, linestyle="--",
                            linewidth=1.0, alpha=0.5)

        damping_note = " (Damped)" if is_damped else " (Undamped)"
        ax.set_title(f"Free Vibration Response{damping_note}")
        ax.set_xlabel("Time  t (s)")
        ax.set_ylabel("Displacement  q(t)")
        ax.legend(loc="upper right", fontsize=9)
        ax.axhline(0, color="black", linewidth=0.8)
        ax.grid(True, alpha=0.35)

        fig.tight_layout()

        if save:
            self._save_fig(fig, "time_response")

        return fig

    # ------------------------------------------------------------------
    # 5. FRF Magnitude Sketch (Asymptotic Bode)
    # ------------------------------------------------------------------

    def plot_frf_magnitude_sketch(
        self,
        frf_result: SolverResult,
        save: bool = False,
    ) -> "plt.Figure":
        """점근선 Bode magnitude 스케치.

        실제 곡선에 직선 점근선을 겹쳐 그려서
        시험의 '스케치' 문제에 대응한다.

        Parameters
        ----------
        frf_result : SolverResult
            FRFSolver 또는 FRFDecomposeSolver 의 출력.
        save : bool
        """
        _apply_style()
        fa = frf_result.final_answer
        given = frf_result.given

        # ---- Build actual FRF ----
        omega, mag_dB = self._build_mag_dB(fa, given)

        # ---- Corner frequencies ----
        cp = fa.get("characteristic_points", {})
        corner_freqs: list[tuple[float, str, int]] = []   # (freq, label, dB/dec slope change)

        for s1 in cp.get("first_order_subsystems", []):
            wb = _safe_float(s1.get("omega_b"))
            if wb and wb > 0:
                corner_freqs.append((wb, f"\u03c9_b={wb:.3g}", -20))

        for s2 in cp.get("second_order_subsystems", []):
            wn = _safe_float(s2.get("omega_n"))
            if wn and wn > 0:
                corner_freqs.append((wn, f"\u03c9_n={wn:.3g}", -40))

        corner_freqs.sort(key=lambda x: x[0])

        # ---- Asymptote construction ----
        dc_gain = _safe_float(fa.get("dc_gain"))
        dc_dB = 20 * np.log10(dc_gain) if dc_gain and dc_gain > 0 else 0.0

        # Build piecewise asymptote over log-frequency grid
        asym_dB = np.full_like(omega, dc_dB)
        for idx, (wc, _, slope_change) in enumerate(corner_freqs):
            mask = omega > wc
            for jdx in range(idx, len(corner_freqs)):
                wc2 = corner_freqs[jdx][0]
                seg_mask = (omega > wc2) if jdx > idx else mask
                delta = slope_change * np.log10(omega[mask] / wc)
                asym_dB[mask] += delta
                # Only apply once per region
                break

        # Recompute properly: cumulative slope changes
        asym_dB2 = np.full_like(omega, dc_dB)
        cumulative_slope = 0
        prev_wc = omega[0]
        sorted_corners = sorted(corner_freqs, key=lambda x: x[0])

        for i, freq_pt in enumerate(omega):
            slope_here = cumulative_slope
            for wc, _, ds in sorted_corners:
                if freq_pt > wc:
                    slope_here += ds
            asym_dB2[i] = dc_dB + slope_here * (np.log10(freq_pt) - np.log10(omega[0]))

        # Reset at first corner
        if sorted_corners:
            w0 = sorted_corners[0][0]
            asym_dB2 = np.full_like(omega, dc_dB)
            for i, freq_pt in enumerate(omega):
                delta = 0.0
                for wc, _, ds in sorted_corners:
                    if freq_pt > wc:
                        delta += ds * (np.log10(freq_pt) - np.log10(wc))
                asym_dB2[i] = dc_dB + delta

        # ---- Figure ----
        fig, ax = plt.subplots(figsize=_FIG_SINGLE)
        fig.suptitle("FRF Magnitude Sketch — Asymptotic Bode Approximation", fontsize=12)

        ax.semilogx(omega, mag_dB, color=_STABLE_COLOR, linewidth=2.0,
                    label="Actual |G(j\u03c9)|", zorder=3)
        ax.semilogx(omega, asym_dB2, color=_UNSTABLE_COLOR, linewidth=1.5,
                    linestyle="--", label="Asymptote", zorder=2)

        # DC gain reference
        ax.axhline(dc_dB, color="gray", linestyle=":", linewidth=1.0,
                   label=f"DC gain = {dc_dB:.1f} dB", alpha=0.7)

        # Annotate corner frequencies
        for wc, label, ds in sorted_corners:
            # Find mag_dB at corner
            idx = np.argmin(np.abs(omega - wc))
            y_val = mag_dB[idx]
            ax.axvline(wc, color=_ACCENT_COLOR, linestyle=":", linewidth=1.2, alpha=0.7)
            ax.annotate(label, xy=(wc, y_val),
                        xytext=(wc * 1.3, y_val + 3),
                        fontsize=8, color=_ACCENT_COLOR,
                        arrowprops=dict(arrowstyle="->", color=_ACCENT_COLOR,
                                        lw=0.8))

        rolloff = _safe_float(
            cp.get("high_freq_rolloff_dB_per_decade")
        )
        if rolloff is not None:
            ax.text(0.97, 0.05, f"High-freq rolloff: {rolloff:.0f} dB/decade",
                    transform=ax.transAxes, ha="right", fontsize=9,
                    color="gray",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                              edgecolor="gray", alpha=0.7))

        ax.set_xlabel("\u03c9 (rad/s)")
        ax.set_ylabel("Magnitude (dB)")
        ax.legend(loc="lower left", fontsize=9)
        ax.grid(True, which="both", alpha=0.3)

        fig.tight_layout()

        if save:
            self._save_fig(fig, "frf_magnitude_sketch")

        return fig

    def _build_mag_dB(
        self, fa: dict, given: dict
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return (omega, mag_dB) arrays, computing from coefficients if needed."""
        omega = fa.get("omega")
        mag_dB = fa.get("magnitude_dB")

        if omega is not None and mag_dB is not None:
            return np.asarray(omega), np.asarray(mag_dB)

        den = given.get("denominator_coefficients") or given.get("coefficients")
        num = given.get("numerator_coefficients") or given.get("numerator", [1])
        if den is None:
            raise ValueError("계수 정보가 없어 Bode 데이터를 계산할 수 없습니다.")

        den_f = [float(c) for c in den]
        num_f = [float(c) for c in num]

        cp = fa.get("characteristic_points", {})
        w_pts = []
        for s1 in cp.get("first_order_subsystems", []):
            if "omega_b" in s1:
                w_pts.append(float(s1["omega_b"]))
        for s2 in cp.get("second_order_subsystems", []):
            if "omega_n" in s2:
                w_pts.append(float(s2["omega_n"]))
        w_min = min(w_pts) / 100.0 if w_pts else 0.01
        w_max = max(w_pts) * 100.0 if w_pts else 1000.0

        omega_arr = np.logspace(np.log10(w_min), np.log10(w_max), 600)
        jw = 1j * omega_arr
        g = np.polyval([complex(c) for c in num_f], jw) / np.polyval([complex(c) for c in den_f], jw)
        return omega_arr, 20 * np.log10(np.abs(g))

    # ------------------------------------------------------------------
    # 6. Stability Diagram
    # ------------------------------------------------------------------

    def plot_stability_diagram(
        self,
        result: SolverResult,
        save: bool = False,
    ) -> "plt.Figure":
        """안정성 다이어그램 — 고유값 위치와 안정/불안정 영역 표시.

        Parameters
        ----------
        result : SolverResult
            StabilitySolver 또는 eigenvalue 정보가 있는 모든 SolverResult.
        save : bool
        """
        _apply_style()
        fa = result.final_answer

        # ---- Collect eigenvalues ----
        raw_eigs = (
            fa.get("eigenvalues_numeric")
            or fa.get("roots")
            or fa.get("eigenvalues")
        )
        if raw_eigs is None:
            raise ValueError(
                "plot_stability_diagram: final_answer에 "
                "'eigenvalues_numeric', 'roots', 'eigenvalues' 키가 없습니다."
            )

        eigs: list[complex] = []
        for e in raw_eigs:
            if isinstance(e, complex):
                eigs.append(e)
            else:
                c = _parse_complex(str(e))
                if c is not None:
                    eigs.append(c)

        if not eigs:
            raise ValueError("plot_stability_diagram: 유효한 고유값 데이터가 없습니다.")

        # ---- Determine stability verdict ----
        stability = fa.get("stability", "Unknown")
        classification = fa.get("classification", "")

        # ---- Figure ----
        fig, ax = plt.subplots(figsize=_FIG_SINGLE)

        # Limits
        margin = max(
            max(abs(e.real) for e in eigs) * 2.5,
            max(abs(e.imag) for e in eigs) * 2.5,
            2.0,
        )

        # ---- Shading ----
        ax.axvspan(-margin, 0, alpha=0.10, color=_NEUTRAL_COLOR,
                   label="Stable region (Re < 0)", zorder=0)
        ax.axvspan(0, margin, alpha=0.10, color=_UNSTABLE_COLOR,
                   label="Unstable region (Re > 0)", zorder=0)

        # ---- Axes ----
        ax.axhline(0, color="black", linewidth=0.8)
        ax.axvline(0, color="black", linewidth=1.2)

        # ---- Plot each eigenvalue ----
        for i, e in enumerate(eigs):
            stable = e.real <= 1e-8
            color = _STABLE_COLOR if stable else _UNSTABLE_COLOR
            ax.plot(e.real, e.imag, "x", color=color,
                    markersize=14, markeredgewidth=3, zorder=5)
            lbl = (
                f"\u03bb\u2081{i+1} = {e.real:.3g}"
                + (f" \u00b1 {abs(e.imag):.3g}j" if abs(e.imag) > 1e-6 else "")
            )
            ax.annotate(
                f" {lbl}",
                xy=(e.real, e.imag),
                fontsize=8, color=color,
                xytext=(e.real + margin * 0.05, e.imag + margin * 0.05),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.5),
            )

        # ---- Verdict box ----
        verdict_color = _NEUTRAL_COLOR if "stable" in stability.lower() else _UNSTABLE_COLOR
        ax.text(0.03, 0.97,
                f"Stability: {stability}\n({classification})",
                transform=ax.transAxes, fontsize=11, va="top",
                color="white",
                bbox=dict(boxstyle="round,pad=0.4", facecolor=verdict_color,
                          edgecolor="none", alpha=0.85))

        ax.set_xlim(-margin, margin)
        ax.set_ylim(-margin, margin)
        ax.set_xlabel("Real part  \u03c3")
        ax.set_ylabel("Imaginary part  j\u03c9")
        ax.set_title(
            f"Stability Diagram — {result.problem_type}", fontsize=12
        )
        ax.legend(loc="lower right", fontsize=9)
        ax.grid(True, alpha=0.30)
        ax.set_aspect("equal", adjustable="datalim")

        fig.tight_layout()

        if save:
            self._save_fig(fig, "stability_diagram")

        return fig

    # ------------------------------------------------------------------
    # Convenience: auto-select plot based on problem_type
    # ------------------------------------------------------------------

    def auto_plot(
        self,
        result: SolverResult,
        save: bool = False,
    ) -> list["plt.Figure"]:
        """problem_type에 따라 적절한 플롯을 자동으로 선택·생성한다.

        Returns
        -------
        list of Figure
            생성된 모든 Figure 객체의 리스트.
        """
        ptype = result.problem_type.lower()
        figs: list["plt.Figure"] = []

        if "frf" in ptype or "frequency response" in ptype or "transfer" in ptype:
            try:
                figs.append(self.plot_bode(result, save=save))
            except Exception as e:
                print(f"  [Bode 플롯 오류] {e}")
            try:
                figs.append(self.plot_frf_magnitude_sketch(result, save=save))
            except Exception as e:
                print(f"  [FRF Sketch 오류] {e}")
            # Pole-zero map
            if "roots" in result.final_answer:
                try:
                    figs.append(self.plot_pole_zero(result, save=save))
                except Exception as e:
                    print(f"  [Pole-Zero 오류] {e}")

        elif "modal" in ptype and "damped" not in ptype:
            try:
                figs.append(self.plot_mode_shapes(result, save=save))
            except Exception as e:
                print(f"  [Mode Shape 오류] {e}")
            try:
                figs.append(self.plot_time_response(result, save=save))
            except Exception as e:
                print(f"  [Time Response 오류] {e}")

        elif "damped" in ptype:
            try:
                figs.append(self.plot_mode_shapes(result, save=save))
            except Exception as e:
                print(f"  [Mode Shape 오류] {e}")
            try:
                figs.append(self.plot_time_response(result, save=save))
            except Exception as e:
                print(f"  [Time Response 오류] {e}")

        elif "stability" in ptype:
            try:
                figs.append(self.plot_stability_diagram(result, save=save))
            except Exception as e:
                print(f"  [Stability 오류] {e}")
            if ("eigenvalues" in result.final_answer
                    or "roots" in result.final_answer):
                try:
                    figs.append(self.plot_pole_zero(result, save=save))
                except Exception as e:
                    print(f"  [Pole-Zero 오류] {e}")

        else:
            # Fallback: try pole-zero if eigenvalues present
            if ("eigenvalues" in result.final_answer
                    or "roots" in result.final_answer):
                try:
                    figs.append(self.plot_pole_zero(result, save=save))
                except Exception as e:
                    print(f"  [Pole-Zero 오류] {e}")

        return figs
