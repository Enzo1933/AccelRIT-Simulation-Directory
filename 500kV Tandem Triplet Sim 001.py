"""
Quadrupole Triplet (FDF) Design Tool for 500kV Tandem Ion Accelerator
======================================================================
Interactive matplotlib tool with sliders for all geometry/beam parameters.
Uses transfer matrix beam optics (hard-edge quad, no fringe fields).

Configuration: F-D-F (focusing-defocusing-focusing in x-plane)
  Q1 (+g1) -> gap1 -> Q2 (-g2) -> gap2 -> Q3 (+g1) -> drift -> target

With two gradient degrees of freedom (g1, g2), the triplet can produce
a round focal spot -- unlike a doublet which is inherently asymmetric.

Dependencies: numpy, matplotlib

Author: Generated for 500kV Tandem Ion Accelerator project at RIT
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox, CheckButtons
from matplotlib.ticker import AutoMinorLocator

# ============================================================
# Constants  (natural / accelerator-physics units)
# ============================================================
PROTON_REST_KEV = 938_272.046     # proton rest energy  [keV]  (= 938.272 MeV)
MU0             = 4 * np.pi * 1e-7  # vacuum permeability [H/m]  (needed for A·t)
C               = 299_792_458.0   # speed of light      [m/s]
IN_TO_M         = 0.0254
MM_TO_M         = 1e-3

# ============================================================
# Physics
# ============================================================
def beam_rigidity(energy_keV):
    """Magnetic rigidity Bρ [T·m] for a proton.

    Natural-units derivation:
      p [keV/c] = sqrt((T + mc²)² - (mc²)²)
      Bρ [T·m]  = p [eV/c] / c [m/s]   (for a singly charged particle)
    No SI mass or charge constants needed.
    """
    p_keV_c = np.sqrt((energy_keV + PROTON_REST_KEV)**2 - PROTON_REST_KEV**2)
    return p_keV_c * 1e3 / C   # keV/c → eV/c, then Bρ = p[eV/c] / c


def quad_transfer_matrix(g, L, Brho):
    """
    Transfer matrix for a magnetic quadrupole (hard-edge model).
    g: field gradient [T/m] (positive = focusing in x, defocusing in y)
    L: effective length [m]
    Brho: magnetic rigidity [T*m]
    Returns (Mx, My) each 2x2 numpy arrays.
    """
    k2 = g / Brho
    if abs(k2) < 1e-12:
        M = np.array([[1, L], [0, 1]])
        return M, M.copy()

    if k2 > 0:
        k = np.sqrt(k2)
        kL = k * L
        Mx = np.array([[np.cos(kL),        np.sin(kL) / k],
                        [-k * np.sin(kL),   np.cos(kL)]])
        My = np.array([[np.cosh(kL),        np.sinh(kL) / k],
                        [k * np.sinh(kL),   np.cosh(kL)]])
    else:
        k = np.sqrt(-k2)
        kL = k * L
        Mx = np.array([[np.cosh(kL),        np.sinh(kL) / k],
                        [k * np.sinh(kL),   np.cosh(kL)]])
        My = np.array([[np.cos(kL),        np.sin(kL) / k],
                        [-k * np.sin(kL),   np.cos(kL)]])
    return Mx, My


def drift_matrix(L):
    M = np.array([[1, L], [0, 1]])
    return M, M.copy()


def total_transfer_matrices(L_outer_m, L_center_m, gap1_m, gap2_m,
                             drift_m, g1, g2, energy_keV):
    """Compute total x and y transfer matrices through the full FDF triplet."""
    Brho = beam_rigidity(energy_keV)
    regions = [
        ('quad',  g1,  L_outer_m),
        ('drift', 0.0, gap1_m),
        ('quad', -g2,  L_center_m),
        ('drift', 0.0, gap2_m),
        ('quad',  g1,  L_outer_m),
        ('drift', 0.0, drift_m),
    ]
    Mx_total = np.eye(2)
    My_total = np.eye(2)
    for rtype, g, length in regions:
        if rtype == 'quad':
            Mx, My = quad_transfer_matrix(g, length, Brho)
        else:
            Mx, My = drift_matrix(length)
        Mx_total = Mx @ Mx_total
        My_total = My @ My_total
    return Mx_total, My_total


def track_envelope(bore_m, L_outer_m, L_center_m, gap1_m, gap2_m,
                   drift_m, g1, g2, energy_keV, x0, xp0, n_steps=400):
    """
    Track beam envelope through FDF triplet.
    Q1(+g1) -> gap1 -> Q2(-g2) -> gap2 -> Q3(+g1) -> drift -> target
    """
    Brho = beam_rigidity(energy_keV)
    total_length = L_outer_m + gap1_m + L_center_m + gap2_m + L_outer_m + drift_m

    regions = [
        ('quad',  g1,  L_outer_m),    # Q1: focusing
        ('drift', 0.0, gap1_m),
        ('quad', -g2,  L_center_m),   # Q2: defocusing
        ('drift', 0.0, gap2_m),
        ('quad',  g1,  L_outer_m),    # Q3: focusing
        ('drift', 0.0, drift_m),
    ]

    z_all = [0.0]
    x_all = [x0]
    y_all = [x0]  # symmetric initial beam
    state_x = np.array([x0, xp0])
    state_y = np.array([x0, xp0])

    for rtype, g, length in regions:
        n = max(int(n_steps * length / total_length), 4)
        dz = length / n
        for _ in range(n):
            if rtype == 'quad':
                Mx, My = quad_transfer_matrix(g, dz, Brho)
            else:
                Mx, My = drift_matrix(dz)
            state_x = Mx @ state_x
            state_y = My @ state_y
            z_all.append(z_all[-1] + dz)
            x_all.append(state_x[0])
            y_all.append(state_y[0])

    z = np.array(z_all)
    x = np.array(x_all)
    y = np.array(y_all)

    def find_crossovers(arr, z_arr):
        crossovers = []
        for i in range(1, len(arr)):
            if arr[i-1] * arr[i] < 0:
                frac = abs(arr[i-1]) / (abs(arr[i-1]) + abs(arr[i]))
                z_cross = z_arr[i-1] + frac * (z_arr[i] - z_arr[i-1])
                crossovers.append(z_cross)
        return crossovers

    q1_end = L_outer_m
    q2_start = L_outer_m + gap1_m
    q2_end = q2_start + L_center_m
    q3_start = q2_end + gap2_m
    q3_end = q3_start + L_outer_m

    return {
        'z': z, 'x': x, 'y': y,
        'final_x': abs(state_x[0]),
        'final_y': abs(state_y[0]),
        'total_length': total_length,
        'q1_end': q1_end,
        'q2_start': q2_start,
        'q2_end': q2_end,
        'q3_start': q3_start,
        'q3_end': q3_end,
        'crossovers_x': find_crossovers(x, z),
        'crossovers_y': find_crossovers(y, z),
        'max_env_x': np.max(np.abs(x)),
        'max_env_y': np.max(np.abs(y)),
    }


def _sweep_gradients(bore_m, L_outer_m, L_center_m, gap1_m, gap2_m,
                     drift_m, energy_keV, x0, xp0, target_r_m,
                     g1_range, g2_range, n_steps, imbalance_thresh):
    """
    Sweep (g1, g2) grid.
    Among candidates with imbalance < thresh, picks closest to target spot.
    Otherwise returns the pair with minimum imbalance.
    """
    best_g1, best_g2 = 0.0, 0.0
    best_imbalance = np.inf
    candidates = []

    for g1 in g1_range:
        for g2 in g2_range:
            env = track_envelope(bore_m, L_outer_m, L_center_m, gap1_m, gap2_m,
                                 drift_m, g1, g2, energy_keV, x0, xp0, n_steps=n_steps)
            if env['max_env_x'] > bore_m * 0.95 or env['max_env_y'] > bore_m * 0.95:
                continue
            fx, fy = env['final_x'], env['final_y']
            avg = (fx + fy) / 2
            if avg < 1e-12:
                continue
            imbalance = abs(fx - fy) / avg

            if imbalance < imbalance_thresh:
                candidates.append((g1, g2, avg, imbalance))
            if imbalance < best_imbalance:
                best_imbalance = imbalance
                best_g1, best_g2 = g1, g2

    if candidates:
        best = min(candidates, key=lambda c: abs(c[2] - target_r_m))
        return best[0], best[1]
    return best_g1, best_g2


def optimize_gradient(bore_m, L_outer_m, L_center_m, gap1_m, gap2_m,
                      drift_m, energy_keV, x0, xp0, target_r_m):
    """
    2D optimization over (g1, g2) for FDF triplet.
    Finds the roundest spot closest to target_r_m.
    """
    # Coarse sweep
    g1, g2 = _sweep_gradients(
        bore_m, L_outer_m, L_center_m, gap1_m, gap2_m, drift_m,
        energy_keV, x0, xp0, target_r_m,
        np.arange(0.5, 25, 1.0), np.arange(0.5, 25, 1.0),
        n_steps=80, imbalance_thresh=0.15)

    if g1 < 0.01:
        return 0.0, 0.0

    # Fine sweep around best
    g1, g2 = _sweep_gradients(
        bore_m, L_outer_m, L_center_m, gap1_m, gap2_m, drift_m,
        energy_keV, x0, xp0, target_r_m,
        np.arange(max(0.1, g1 - 1.5), g1 + 1.5, 0.1),
        np.arange(max(0.1, g2 - 1.5), g2 + 1.5, 0.1),
        n_steps=200, imbalance_thresh=0.10)

    return g1, g2


# ============================================================
# Unit conversion helpers
# ============================================================
class UnitState:
    """Manages display units (mm or inches). Internal state always in inches."""
    def __init__(self):
        self.use_inches = True

    @property
    def label(self):
        return "in" if self.use_inches else "mm"

    def to_display(self, val_in):
        return val_in if self.use_inches else val_in * 25.4

    def from_display(self, val_disp):
        return val_disp if self.use_inches else val_disp / 25.4

    def m_to_display(self, val_m):
        return val_m / IN_TO_M if self.use_inches else val_m / MM_TO_M

    def display_to_m(self, val_disp):
        return val_disp * IN_TO_M if self.use_inches else val_disp * MM_TO_M


# ============================================================
# Pre-generated particle distribution (uniform, zero-emittance, correlated)
# Angle proportional to position (matches single-ray envelope model).
# Fixed seed so heatmap doesn't flicker on slider updates.
# ============================================================
N_PARTICLES = 200000
_rng = np.random.default_rng(42)
_r = np.sqrt(_rng.uniform(0, 1, N_PARTICLES))  # radial (uniform in area)
_theta = _rng.uniform(0, 2 * np.pi, N_PARTICLES)
UNIT_X  = _r * np.cos(_theta)   # multiply by x0 -> x position
UNIT_Y  = _r * np.sin(_theta)   # multiply by x0 -> y position
UNIT_XP = UNIT_X                 # multiply by xp0 -> x angle (proportional to x)
UNIT_YP = UNIT_Y                 # multiply by xp0 -> y angle (proportional to y)


# ============================================================
# Main interactive figure
# ============================================================
def main():
    units = UnitState()

    # --- Default parameters ---
    params = {
        'bore_in':      3.0,
        'l_outer_in':   3.0,   # Q1, Q3 length [inches]
        'l_center_in':  4.0,   # Q2 length [inches]
        'gap1_in':      2.0,   # Q1-Q2 edge-to-edge gap [inches]
        'gap2_in':      2.0,   # Q2-Q3 edge-to-edge gap [inches]
        'drift_in':    18.0,   # drift to target [inches]
        'target_mm':    2.0,   # target spot radius [mm]
        'energy_keV':1000.0,   # beam kinetic energy [keV]
        'x0_in':        0.20,  # initial beam radius [inches]
        'xp0_mrad':    30.0,   # initial divergence [mrad]
        'auto_opt':     True,
        'manual_g1':    5.0,   # manual outer gradient [T/m]
        'manual_g2':   10.0,   # manual center gradient [T/m]
    }

    # ---- Figure layout ----
    fig = plt.figure(figsize=(14, 10))
    fig.patch.set_facecolor('#f8f8f6')
    fig.suptitle('Quadrupole Triplet (FDF) Design Tool — 500kV Tandem Ion Accelerator',
                 fontsize=14, fontweight='bold', y=0.98)

    ax = fig.add_axes([0.06, 0.44, 0.52, 0.48])
    ax_heatmap = fig.add_axes([0.66, 0.46, 0.24, 0.44])
    ax_cbar = fig.add_axes([0.91, 0.46, 0.015, 0.44])

    # Slider layout constants
    slider_left   = 0.12
    slider_width  = 0.22
    slider_height = 0.018
    slider_x2     = 0.58
    tb_width      = 0.055
    tb_gap        = 0.008
    row_step      = 0.035
    row0_y        = 0.36

    def make_slider_ax(x, row):
        y = row0_y - row * row_step
        return fig.add_axes([x, y, slider_width, slider_height])

    def make_tb_ax(x, row):
        y = row0_y - row * row_step
        return fig.add_axes([x + slider_width + tb_gap, y - 0.001, tb_width, slider_height + 0.002])

    # ---- Column 1: geometry sliders (rows 0-7) ----
    s_bore     = Slider(make_slider_ax(slider_left, 0), 'Bore radius (in)',        0.25,  6.0, valinit=params['bore_in'],     valstep=0.1)
    s_l_outer  = Slider(make_slider_ax(slider_left, 1), 'Outer quad L (in)',       1.0,   6.0, valinit=params['l_outer_in'],  valstep=0.1)
    s_l_center = Slider(make_slider_ax(slider_left, 2), 'Center quad L (in)',      1.0,   8.0, valinit=params['l_center_in'], valstep=0.1)
    s_gap1     = Slider(make_slider_ax(slider_left, 3), 'Q1-Q2 gap (in)',          0.5,   8.0, valinit=params['gap1_in'],     valstep=0.1)
    s_gap2     = Slider(make_slider_ax(slider_left, 4), 'Q2-Q3 gap (in)',          0.5,   8.0, valinit=params['gap2_in'],     valstep=0.1)
    s_drift    = Slider(make_slider_ax(slider_left, 5), 'Drift to target (in)',    6.0,  36.0, valinit=params['drift_in'],    valstep=0.5)
    s_x0       = Slider(make_slider_ax(slider_left, 6), 'Init spot r (in)',        0.01,  0.25, valinit=params['x0_in'],      valstep=0.01)
    s_xp0      = Slider(make_slider_ax(slider_left, 7), 'Init divergence (mrad)',  0.5, 100.0, valinit=params['xp0_mrad'],    valstep=0.5)

    # ---- Column 2: beam/gradient sliders (rows 0-3) ----
    s_energy = Slider(make_slider_ax(slider_x2, 0), 'Beam energy (keV)',  100.0, 3000.0, valinit=params['energy_keV'], valstep=50.0)
    s_spot   = Slider(make_slider_ax(slider_x2, 1), 'Target spot (mm)',     0.1,   10.0, valinit=params['target_mm'],  valstep=0.1)
    s_g1     = Slider(make_slider_ax(slider_x2, 2), 'Grad outer (T/m)',     0.1,   30.0, valinit=params['manual_g1'],  valstep=0.1)
    s_g2     = Slider(make_slider_ax(slider_x2, 3), 'Grad center (T/m)',    0.1,   30.0, valinit=params['manual_g2'],  valstep=0.1)

    all_sliders = [s_bore, s_l_outer, s_l_center, s_gap1, s_gap2, s_drift, s_x0, s_xp0,
                   s_energy, s_spot, s_g1, s_g2]
    for s in all_sliders:
        s.label.set_fontsize(9)
        s.valtext.set_fontsize(9)

    # ---- Text input boxes ----
    tb_bore     = TextBox(make_tb_ax(slider_left, 0), '', initial=f'{params["bore_in"]:.1f}')
    tb_l_outer  = TextBox(make_tb_ax(slider_left, 1), '', initial=f'{params["l_outer_in"]:.1f}')
    tb_l_center = TextBox(make_tb_ax(slider_left, 2), '', initial=f'{params["l_center_in"]:.1f}')
    tb_gap1     = TextBox(make_tb_ax(slider_left, 3), '', initial=f'{params["gap1_in"]:.1f}')
    tb_gap2     = TextBox(make_tb_ax(slider_left, 4), '', initial=f'{params["gap2_in"]:.1f}')
    tb_drift    = TextBox(make_tb_ax(slider_left, 5), '', initial=f'{params["drift_in"]:.1f}')
    tb_x0       = TextBox(make_tb_ax(slider_left, 6), '', initial=f'{params["x0_in"]:.2f}')
    tb_xp0      = TextBox(make_tb_ax(slider_left, 7), '', initial=f'{params["xp0_mrad"]:.1f}')
    tb_energy   = TextBox(make_tb_ax(slider_x2,   0), '', initial=f'{params["energy_keV"]:.0f}')
    tb_spot     = TextBox(make_tb_ax(slider_x2,   1), '', initial=f'{params["target_mm"]:.1f}')
    tb_g1       = TextBox(make_tb_ax(slider_x2,   2), '', initial=f'{params["manual_g1"]:.2f}')
    tb_g2       = TextBox(make_tb_ax(slider_x2,   3), '', initial=f'{params["manual_g2"]:.2f}')

    def connect_slider_textbox(slider, textbox, fmt):
        def on_submit(text):
            try:
                val = float(np.clip(float(text), slider.valmin, slider.valmax))
                slider.set_val(val)
            except ValueError:
                pass
            textbox.set_val(fmt.format(slider.val))

        def on_slider_changed(val):
            textbox.set_val(fmt.format(val))

        textbox.on_submit(on_submit)
        slider.on_changed(on_slider_changed)

    connect_slider_textbox(s_bore,     tb_bore,     '{:.1f}')
    connect_slider_textbox(s_l_outer,  tb_l_outer,  '{:.1f}')
    connect_slider_textbox(s_l_center, tb_l_center, '{:.1f}')
    connect_slider_textbox(s_gap1,     tb_gap1,     '{:.1f}')
    connect_slider_textbox(s_gap2,     tb_gap2,     '{:.1f}')
    connect_slider_textbox(s_drift,    tb_drift,    '{:.1f}')
    connect_slider_textbox(s_x0,       tb_x0,       '{:.2f}')
    connect_slider_textbox(s_xp0,      tb_xp0,      '{:.1f}')
    connect_slider_textbox(s_energy,   tb_energy,   '{:.0f}')
    connect_slider_textbox(s_spot,     tb_spot,     '{:.1f}')
    connect_slider_textbox(s_g1,       tb_g1,       '{:.2f}')
    connect_slider_textbox(s_g2,       tb_g2,       '{:.2f}')

    # ---- Gray out gradient sliders in auto mode ----
    def set_grad_active(active):
        alpha = 1.0 if active else 0.3
        for s, tb in [(s_g1, tb_g1), (s_g2, tb_g2)]:
            s.ax.set_alpha(alpha)
            s.active = active
            tb.ax.set_alpha(alpha)
        fig.canvas.draw_idle()

    # ---- Auto-optimize checkbox (col 2, row 4) ----
    ax_check = fig.add_axes([slider_x2, row0_y - 4 * row_step, 0.2, 0.025])
    check_auto = CheckButtons(ax_check, ['Auto-optimize gradients'], [params['auto_opt']])
    check_auto.labels[0].set_fontsize(9)
    set_grad_active(not params['auto_opt'])

    # ---- Unit toggle button ----
    ax_unit_btn = fig.add_axes([slider_x2 + 0.22, row0_y - 4 * row_step, 0.06, 0.025])
    btn_unit = Button(ax_unit_btn, 'mm <> in', hovercolor='#ddd')
    btn_unit.label.set_fontsize(8)

    # ---- Info text area ----
    ax_info = fig.add_axes([0.58, 0.01, 0.40, 0.17])
    ax_info.set_axis_off()
    info_text = ax_info.text(0, 1, '', fontsize=8, fontfamily='monospace',
                              verticalalignment='top', transform=ax_info.transAxes)

    # ---- Beam info (bottom left) ----
    ax_beam = fig.add_axes([0.08, 0.02, 0.45, 0.06])
    ax_beam.set_axis_off()
    beam_text = ax_beam.text(0, 1, '', fontsize=8.5, fontfamily='monospace',
                              verticalalignment='top', transform=ax_beam.transAxes,
                              color='#666')

    # ---- Update function ----
    def update(val=None):
        ax.clear()

        bore_in     = s_bore.val
        l_outer_in  = s_l_outer.val
        l_center_in = s_l_center.val
        gap1_in     = s_gap1.val
        gap2_in     = s_gap2.val
        drift_in    = s_drift.val
        x0_in       = s_x0.val
        xp0_mrad    = s_xp0.val
        target_mm   = s_spot.val
        energy_keV  = s_energy.val
        manual_g1   = s_g1.val
        manual_g2   = s_g2.val
        auto_opt    = params['auto_opt']

        bore_m     = bore_in * IN_TO_M
        l_outer_m  = l_outer_in * IN_TO_M
        l_center_m = l_center_in * IN_TO_M
        gap1_m     = gap1_in * IN_TO_M
        gap2_m     = gap2_in * IN_TO_M
        drift_m    = drift_in * IN_TO_M
        x0_m       = x0_in * IN_TO_M
        xp0_rad    = xp0_mrad * 1e-3
        target_m   = target_mm * MM_TO_M
        Brho       = beam_rigidity(energy_keV)

        # Get gradients
        if auto_opt:
            g1, g2 = optimize_gradient(bore_m, l_outer_m, l_center_m, gap1_m, gap2_m,
                                        drift_m, energy_keV, x0_m, xp0_rad, target_m)
        else:
            g1, g2 = manual_g1, manual_g2

        if g1 < 0.01:
            ax.text(0.5, 0.5,
                    'No solution found.\n'
                    'Try larger bore, shorter drift, or different gap sizes.',
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=12, color='red')
            info_text.set_text('')
            fig.canvas.draw_idle()
            return

        env = track_envelope(bore_m, l_outer_m, l_center_m, gap1_m, gap2_m,
                              drift_m, g1, g2, energy_keV, x0_m, xp0_rad)

        z = env['z']
        x = env['x']
        y_env = env['y']

        # Convert to display units
        z_d = np.array([units.m_to_display(v) for v in z])
        x_d = np.array([units.m_to_display(v) for v in x])
        y_d = np.array([units.m_to_display(v) for v in y_env])
        bore_d   = units.m_to_display(bore_m)
        target_d = units.m_to_display(target_m)
        du = units.label

        # Region boundaries in display units
        q1_end_d   = units.m_to_display(env['q1_end'])
        q2_start_d = units.m_to_display(env['q2_start'])
        q2_end_d   = units.m_to_display(env['q2_end'])
        q3_start_d = units.m_to_display(env['q3_start'])
        q3_end_d   = units.m_to_display(env['q3_end'])
        total_d    = units.m_to_display(env['total_length'])

        # ---- Plot ----
        # Quad shading (blue for focusing, red-ish for defocusing)
        ax.axvspan(0, q1_end_d, alpha=0.08, color='#4444cc')
        ax.axvspan(q2_start_d, q2_end_d, alpha=0.08, color='#cc4444')
        ax.axvspan(q3_start_d, q3_end_d, alpha=0.08, color='#4444cc')
        ax.text(q1_end_d / 2, bore_d * 0.92, 'Q1 (F)',
                ha='center', fontsize=9, color='#888', style='italic')
        ax.text((q2_start_d + q2_end_d) / 2, bore_d * 0.92, 'Q2 (D)',
                ha='center', fontsize=9, color='#888', style='italic')
        ax.text((q3_start_d + q3_end_d) / 2, bore_d * 0.92, 'Q3 (F)',
                ha='center', fontsize=9, color='#888', style='italic')

        # Bore limits
        ax.axhline(bore_d, color='red', ls='--', lw=0.8, alpha=0.35,
                   label=f'Bore +/-{bore_in:.1f} in')
        ax.axhline(-bore_d, color='red', ls='--', lw=0.8, alpha=0.35)

        # Target spot lines
        ax.axhline(target_d, color='gray', ls=':', lw=0.6, alpha=0.5,
                   label=f'Target +/-{target_mm:.1f} mm')
        ax.axhline(-target_d, color='gray', ls=':', lw=0.6, alpha=0.5)

        # Beam envelopes
        ax.plot(z_d, x_d, color='#2266bb', lw=2, label='x-plane')
        ax.plot(z_d, -x_d, color='#2266bb', lw=2)
        ax.plot(z_d, y_d, color='#cc5522', lw=2, label='y-plane')
        ax.plot(z_d, -y_d, color='#cc5522', lw=2)

        ax.axhline(0, color='black', lw=0.5, alpha=0.3)

        # Crossover markers
        for zc in env['crossovers_x']:
            ax.plot(units.m_to_display(zc), 0, 'o', color='#2266bb', ms=5, zorder=5)
        for zc in env['crossovers_y']:
            ax.plot(units.m_to_display(zc), 0, 's', color='#cc5522', ms=5, zorder=5)

        ax.set_xlabel(f'z ({du})', fontsize=10)
        ax.set_ylabel(f'r ({du})', fontsize=10)
        ax.set_xlim(-total_d * 0.02, total_d * 1.04)
        max_r = max(bore_d, np.max(np.abs(x_d)), np.max(np.abs(y_d))) * 1.15
        ax.set_ylim(-max_r, max_r)
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        ax.tick_params(labelsize=8)
        ax.set_title(f'Beam envelope — {energy_keV:.0f} keV proton  |  '
                     f'g_outer = {g1:.2f}, g_center = {g2:.2f} T/m',
                     fontsize=10, pad=8)
        ax.grid(True, alpha=0.15)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(which='minor', length=2.5, color='gray')
        ax.grid(True, which='minor', alpha=0.06)

        # Clipping warning
        clipping = env['max_env_x'] > bore_m * 0.95 or env['max_env_y'] > bore_m * 0.95
        if clipping:
            ax.text(0.5, 0.95, 'BEAM CLIPPED BY BORE',
                    ha='center', va='top', transform=ax.transAxes,
                    fontsize=11, color='red', fontweight='bold',
                    bbox=dict(boxstyle='round', facecolor='#fee', alpha=0.9))

        # ---- Compute outputs ----
        pole_tip_outer  = g1 * bore_m
        pole_tip_center = g2 * bore_m
        at_outer  = pole_tip_outer * bore_m / MU0
        at_center = pole_tip_center * bore_m / MU0

        final_x_mm = env['final_x'] * 1000
        final_y_mm = env['final_y'] * 1000
        avg_mm = (final_x_mm + final_y_mm) / 2
        imbalance = abs(final_x_mm - final_y_mm) / avg_mm * 100 if avg_mm > 0 else 0

        cx_str = ', '.join([f"{units.m_to_display(zc):.1f}" for zc in env['crossovers_x']])
        cy_str = ', '.join([f"{units.m_to_display(zc):.1f}" for zc in env['crossovers_y']])
        if not cx_str: cx_str = 'none'
        if not cy_str: cy_str = 'none'

        fx_d = units.m_to_display(env['final_x'])
        fy_d = units.m_to_display(env['final_y'])
        avg_d = (fx_d + fy_d) / 2
        spot_dec = 4 if units.use_inches else 2

        info_lines = [
            f"{'RESULTS':─^44}",
            f"  g_outer (Q1,Q3): {g1:>6.2f} T/m   g_center (Q2): {g2:>6.2f} T/m",
            f"  B_tip outer: {pole_tip_outer*1000:>7.1f} mT   center: {pole_tip_center*1000:>7.1f} mT",
            f"  A·t  outer:  {at_outer:>7.0f}       center: {at_center:>7.0f}",
            f"  {'':─^42}",
            f"  Spot radius x:  {fx_d:>8.{spot_dec}f}  {du}  ({final_x_mm:.2f} mm)",
            f"  Spot radius y:  {fy_d:>8.{spot_dec}f}  {du}  ({final_y_mm:.2f} mm)",
            f"  Avg spot radius: {avg_d:>7.{spot_dec}f}  {du}  ({avg_mm:.2f} mm)",
            f"  x/y imbalance:   {imbalance:>7.1f}  %",
            f"  {'':─^42}",
            f"  x crossovers (z): {cx_str} {du}",
            f"  y crossovers (z): {cy_str} {du}",
        ]
        info_text.set_text('\n'.join(info_lines))

        beam_lines = [
            f"Beam: {energy_keV:.0f} keV proton  |  "
            f"Brho = {Brho:.4f} T·m  |  "
            f"r\u2080 = {x0_m*1000:.1f} mm  |  "
            f"\u03b8\u2080 = {xp0_mrad:.1f} mrad",
        ]
        beam_text.set_text('\n'.join(beam_lines))

        # ---- Target cross-section heatmap ----
        ax_heatmap.clear()
        ax_cbar.clear()

        Mx_total, My_total = total_transfer_matrices(
            l_outer_m, l_center_m, gap1_m, gap2_m, drift_m, g1, g2, energy_keV)

        x_final = Mx_total[0, 0] * (x0_m * UNIT_X) + Mx_total[0, 1] * (xp0_rad * UNIT_XP)
        y_final = My_total[0, 0] * (x0_m * UNIT_Y) + My_total[0, 1] * (xp0_rad * UNIT_YP)

        x_mm = x_final * 1000
        y_mm = y_final * 1000

        env_x_mm = env['final_x'] * 1000
        env_y_mm = env['final_y'] * 1000
        avg_spot_mm = (env_x_mm + env_y_mm) / 2

        # Extent based on actual particle spread and target, whichever is larger
        particle_max = max(np.max(np.abs(x_mm)), np.max(np.abs(y_mm)))
        extent_mm = max(particle_max, target_mm) * 1.4
        if extent_mm < 0.1:
            extent_mm = 1.0

        n_bins = 80
        H, _, _ = np.histogram2d(
            x_mm, y_mm, bins=n_bins,
            range=[[-extent_mm, extent_mm], [-extent_mm, extent_mm]])
        H = H.T  # imshow expects (row, col) = (y, x)

        im = ax_heatmap.imshow(
            H, extent=[-extent_mm, extent_mm, -extent_mm, extent_mm],
            origin='lower', cmap='inferno', aspect='equal',
            interpolation='bilinear')
        fig.colorbar(im, cax=ax_cbar, label='Particle density')
        ax_cbar.yaxis.label.set_fontsize(8)
        ax_cbar.tick_params(labelsize=7)

        # Circle overlays
        theta_circ = np.linspace(0, 2 * np.pi, 100)
        ax_heatmap.plot(target_mm * np.cos(theta_circ), target_mm * np.sin(theta_circ),
                        'c--', lw=1, alpha=0.8, label=f'Target {target_mm:.1f} mm')
        ax_heatmap.plot(avg_spot_mm * np.cos(theta_circ), avg_spot_mm * np.sin(theta_circ),
                        'w--', lw=1, alpha=0.8, label=f'Actual {avg_spot_mm:.2f} mm')
        ax_heatmap.legend(loc='upper right', fontsize=6, framealpha=0.6)

        ax_heatmap.set_xlabel('x (mm)', fontsize=9)
        ax_heatmap.set_ylabel('y (mm)', fontsize=9)
        ax_heatmap.set_title('Target cross-section', fontsize=10)
        ax_heatmap.tick_params(labelsize=7)

        fig.canvas.draw_idle()

    # ---- Connect sliders ----
    for s in all_sliders:
        s.on_changed(update)

    def on_check(label):
        params['auto_opt'] = not params['auto_opt']
        set_grad_active(not params['auto_opt'])
        update()
    check_auto.on_clicked(on_check)

    def on_unit_toggle(event):
        units.use_inches = not units.use_inches
        new_label = units.label

        # Only spatial geometry sliders change with the unit toggle
        s_bore.label.set_text(f'Bore radius ({new_label})')
        s_l_outer.label.set_text(f'Outer quad L ({new_label})')
        s_l_center.label.set_text(f'Center quad L ({new_label})')
        s_gap1.label.set_text(f'Q1-Q2 gap ({new_label})')
        s_gap2.label.set_text(f'Q2-Q3 gap ({new_label})')
        s_drift.label.set_text(f'Drift to target ({new_label})')
        s_x0.label.set_text(f'Init spot r ({new_label})')
        # s_xp0 stays in mrad, s_energy stays in keV, s_spot stays in mm

        if units.use_inches:
            s_bore.valmin, s_bore.valmax, s_bore.valstep           = 0.25, 6.0, 0.1
            s_l_outer.valmin, s_l_outer.valmax, s_l_outer.valstep  = 1.0, 6.0, 0.1
            s_l_center.valmin, s_l_center.valmax, s_l_center.valstep = 1.0, 8.0, 0.1
            s_gap1.valmin, s_gap1.valmax, s_gap1.valstep            = 0.5, 8.0, 0.1
            s_gap2.valmin, s_gap2.valmax, s_gap2.valstep            = 0.5, 8.0, 0.1
            s_drift.valmin, s_drift.valmax, s_drift.valstep         = 6.0, 36.0, 0.5
            s_x0.valmin, s_x0.valmax, s_x0.valstep                 = 0.01, 0.25, 0.01
        else:
            s_bore.valmin, s_bore.valmax, s_bore.valstep           = 6.0, 152.0, 1.0
            s_l_outer.valmin, s_l_outer.valmax, s_l_outer.valstep  = 25.0, 152.0, 1.0
            s_l_center.valmin, s_l_center.valmax, s_l_center.valstep = 25.0, 203.0, 1.0
            s_gap1.valmin, s_gap1.valmax, s_gap1.valstep            = 13.0, 203.0, 1.0
            s_gap2.valmin, s_gap2.valmax, s_gap2.valstep            = 13.0, 203.0, 1.0
            s_drift.valmin, s_drift.valmax, s_drift.valstep         = 152.0, 914.0, 5.0
            s_x0.valmin, s_x0.valmax, s_x0.valstep                 = 0.25, 6.35, 0.1

        update()

    btn_unit.on_clicked(on_unit_toggle)

    # ---- Initial draw ----
    set_grad_active(not params['auto_opt'])
    update()

    Brho_init = beam_rigidity(params['energy_keV'])
    print("=" * 60)
    print("  Quadrupole Triplet (FDF) Design Tool")
    print(f"  {params['energy_keV']:.0f} keV Proton | F-D-F configuration")
    print("  Drag sliders or type values. Close window to exit.")
    print("=" * 60)
    print(f"  Brho = {Brho_init:.4f} T·m  at {params['energy_keV']:.0f} keV")
    print(f"  2 DOF (g_outer, g_center) -> can achieve round spot")
    print()

    plt.show()


if __name__ == '__main__':
    main()
