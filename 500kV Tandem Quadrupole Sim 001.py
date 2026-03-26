"""
Quadrupole Doublet Design Tool for 1 MeV Proton Beam Focusing
=============================================================
Interactive matplotlib tool with sliders for all geometry/beam parameters.
Uses transfer matrix beam optics (hard-edge quad, no fringe fields).

Dependencies: numpy, matplotlib
Run: python quad_doublet_tool.py

Author: Generated for 500kV Tandem Ion Accelerator project at RIT
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox, CheckButtons
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import AutoMinorLocator

# ============================================================
# Constants
# ============================================================
PROTON_MASS_KG = 1.6726219e-27
PROTON_CHARGE  = 1.602176634e-19
MU0            = 4 * np.pi * 1e-7
MEV_TO_J       = 1.602176634e-13
C              = 299792458.0
IN_TO_M        = 0.0254
MM_TO_M        = 1e-3

# ============================================================
# Physics
# ============================================================
def beam_rigidity(energy_MeV):
    """Compute magnetic rigidity Bρ [T·m] for a proton at given kinetic energy."""
    E_J = energy_MeV * MEV_TO_J
    E_rest = PROTON_MASS_KG * C**2
    p = np.sqrt((E_J + E_rest)**2 - E_rest**2) / C
    return p / PROTON_CHARGE


def quad_transfer_matrix(g, L, Brho):
    """
    Transfer matrix for a magnetic quadrupole (hard-edge model).
    g: field gradient [T/m] (positive = focusing in x, defocusing in y)
    L: effective length [m]
    Brho: magnetic rigidity [T·m]
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


def track_envelope(bore_m, L_mag_m, gap_m, drift_m, g1, energy_MeV, x0, xp0, n_steps=400):
    """
    Track beam envelope through FD doublet.
    Returns dict with z positions, x/y envelopes, region boundaries, crossovers, etc.
    """
    Brho = beam_rigidity(energy_MeV)
    total_length = L_mag_m + gap_m + L_mag_m + drift_m

    # Define regions: Q1(+g), drift, Q2(-g), drift to target
    regions = [
        ('quad',  g1,  L_mag_m),
        ('drift', 0.0, gap_m),
        ('quad', -g1,  L_mag_m),
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

    # Find zero crossings via sign changes
    def find_crossovers(arr, z_arr):
        crossovers = []
        for i in range(1, len(arr)):
            if arr[i-1] * arr[i] < 0:
                frac = abs(arr[i-1]) / (abs(arr[i-1]) + abs(arr[i]))
                z_cross = z_arr[i-1] + frac * (z_arr[i] - z_arr[i-1])
                crossovers.append(z_cross)
        return crossovers

    return {
        'z': z, 'x': x, 'y': y,
        'final_x': abs(state_x[0]),
        'final_y': abs(state_y[0]),
        'total_length': total_length,
        'q1_end': L_mag_m,
        'q2_start': L_mag_m + gap_m,
        'q2_end': 2 * L_mag_m + gap_m,
        'crossovers_x': find_crossovers(x, z),
        'crossovers_y': find_crossovers(y, z),
        'max_env_x': np.max(np.abs(x)),
        'max_env_y': np.max(np.abs(y)),
    }


def optimize_gradient(bore_m, L_mag_m, gap_m, drift_m, energy_MeV, x0, xp0,
                       target_r_m, balance_tol=0.30):
    """
    Find optimal gradient that:
      1. Produces avg spot size closest to target_r_m
      2. Keeps x and y spot radii within balance_tol (30%) of each other
      3. Doesn't clip the bore
    """
    best_g = 0.0
    best_score = np.inf

    # Coarse sweep
    for g1 in np.arange(0.1, 50, 0.25):
        env = track_envelope(bore_m, L_mag_m, gap_m, drift_m, g1, energy_MeV, x0, xp0, n_steps=200)
        if env['max_env_x'] > bore_m * 0.95 or env['max_env_y'] > bore_m * 0.95:
            continue
        fx, fy = env['final_x'], env['final_y']
        avg = (fx + fy) / 2
        if avg < 1e-12:
            continue

        # Check balance: |fx - fy| / avg < tolerance
        imbalance = abs(fx - fy) / avg
        if imbalance > balance_tol:
            continue

        score = abs(avg - target_r_m)
        if score < best_score:
            best_score = score
            best_g = g1

    # Fine sweep around best
    if best_g > 0:
        for g1 in np.arange(max(0.05, best_g - 0.5), best_g + 0.5, 0.01):
            env = track_envelope(bore_m, L_mag_m, gap_m, drift_m, g1, energy_MeV, x0, xp0, n_steps=200)
            if env['max_env_x'] > bore_m * 0.95 or env['max_env_y'] > bore_m * 0.95:
                continue
            fx, fy = env['final_x'], env['final_y']
            avg = (fx + fy) / 2
            if avg < 1e-12:
                continue
            imbalance = abs(fx - fy) / avg
            if imbalance > balance_tol:
                continue
            score = abs(avg - target_r_m)
            if score < best_score:
                best_score = score
                best_g = g1

    return best_g


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
        """Convert meters to display units."""
        return val_m / IN_TO_M if self.use_inches else val_m / MM_TO_M

    def display_to_m(self, val_disp):
        return val_disp * IN_TO_M if self.use_inches else val_disp * MM_TO_M

    def spot_to_display(self, val_mm):
        """Spot size stored in mm internally."""
        return val_mm / 25.4 if self.use_inches else val_mm

    def display_to_spot(self, val_disp):
        return val_disp * 25.4 if self.use_inches else val_disp


# ============================================================
# Main interactive figure
# ============================================================
def main():
    units = UnitState()

    # --- Default parameters (stored in inches / mm internally) ---
    params = {
        'bore_in':    3.0,   # bore radius [inches]
        'lmag_in':    4.0,   # magnet length [inches]
        'spacing_in': 8.0,   # Q1-Q2 center-to-center [inches]
        'drift_in':   18.0,  # drift to target [inches]
        'target_mm':  5.0,   # target spot radius [mm]
        'energy':     1.0,   # beam energy [MeV]
        'x0':         0.005, # initial beam radius [m]
        'xp0':        0.030, # initial divergence [rad]
        'auto_opt':   True,
        'manual_g':   5.0,   # manual gradient [T/m]
    }

    Brho = beam_rigidity(params['energy'])

    # ---- Figure layout ----
    fig = plt.figure(figsize=(14, 9))
    fig.patch.set_facecolor('#f8f8f6')
    fig.suptitle('Quadrupole Doublet Design Tool — 1 MeV Proton Beam',
                 fontsize=14, fontweight='bold', y=0.98)

    # Main plot axes
    ax = fig.add_axes([0.08, 0.38, 0.88, 0.52])

    # Slider region
    slider_left   = 0.12
    slider_width  = 0.22   # shortened to leave room for text input boxes
    slider_height = 0.018
    slider_x2     = 0.58   # second column
    tb_width      = 0.055  # textbox width
    tb_gap        = 0.008  # gap between slider right edge and textbox

    # ---- Create sliders ----
    def make_slider_ax(x, row):
        y = 0.28 - row * 0.04
        return fig.add_axes([x, y, slider_width, slider_height])

    def make_tb_ax(x, row):
        """Textbox axes positioned just to the right of the corresponding slider."""
        y = 0.28 - row * 0.04
        return fig.add_axes([x + slider_width + tb_gap, y - 0.001, tb_width, slider_height + 0.002])

    ax_bore    = make_slider_ax(slider_left, 0)
    ax_lmag    = make_slider_ax(slider_left, 1)
    ax_spacing = make_slider_ax(slider_left, 2)
    ax_drift   = make_slider_ax(slider_left, 3)
    ax_spot    = make_slider_ax(slider_x2, 0)
    ax_grad    = make_slider_ax(slider_x2, 1)

    s_bore = Slider(ax_bore, 'Bore radius (in)', 1.5, 6.0, valinit=params['bore_in'], valstep=0.1)
    s_lmag = Slider(ax_lmag, 'Mag length (in)', 2.0, 6.0, valinit=params['lmag_in'], valstep=0.1)
    s_spacing = Slider(ax_spacing, 'Q1–Q2 spacing (in)', 3.0, 12.0, valinit=params['spacing_in'], valstep=0.1)
    s_drift = Slider(ax_drift, 'Drift to target (in)', 12.0, 24.0, valinit=params['drift_in'], valstep=0.5)
    s_spot = Slider(ax_spot, 'Target spot (mm)', 0.1, 10.0, valinit=params['target_mm'], valstep=0.1)
    s_grad = Slider(ax_grad, 'Gradient (T/m)', 0.1, 40.0, valinit=params['manual_g'], valstep=0.1)

    all_sliders = [s_bore, s_lmag, s_spacing, s_drift, s_spot, s_grad]
    for s in all_sliders:
        s.label.set_fontsize(9)
        s.valtext.set_fontsize(9)

    # ---- Text input boxes (one per slider, right of each slider) ----
    # User can type an exact value and press Enter; slider and plot update immediately.
    tb_bore    = TextBox(make_tb_ax(slider_left, 0), '', initial=f'{params["bore_in"]:.1f}')
    tb_lmag    = TextBox(make_tb_ax(slider_left, 1), '', initial=f'{params["lmag_in"]:.1f}')
    tb_spacing = TextBox(make_tb_ax(slider_left, 2), '', initial=f'{params["spacing_in"]:.1f}')
    tb_drift   = TextBox(make_tb_ax(slider_left, 3), '', initial=f'{params["drift_in"]:.1f}')
    tb_spot    = TextBox(make_tb_ax(slider_x2, 0),   '', initial=f'{params["target_mm"]:.1f}')
    tb_grad    = TextBox(make_tb_ax(slider_x2, 1),   '', initial=f'{params["manual_g"]:.2f}')

    def connect_slider_textbox(slider, textbox, fmt):
        """
        Wire a slider and textbox bidirectionally.
          - Enter in textbox  → clamp to slider range → set_val → triggers update()
          - Slider drag       → update textbox display (set_val does NOT fire on_submit)
        """
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

    connect_slider_textbox(s_bore,    tb_bore,    '{:.1f}')
    connect_slider_textbox(s_lmag,    tb_lmag,    '{:.1f}')
    connect_slider_textbox(s_spacing, tb_spacing, '{:.1f}')
    connect_slider_textbox(s_drift,   tb_drift,   '{:.1f}')
    connect_slider_textbox(s_spot,    tb_spot,    '{:.1f}')
    connect_slider_textbox(s_grad,    tb_grad,    '{:.2f}')

    # Gray out gradient slider AND its textbox when auto-optimize is on
    def set_grad_active(active):
        alpha = 1.0 if active else 0.3
        s_grad.ax.set_alpha(alpha)
        s_grad.active = active
        tb_grad.ax.set_alpha(alpha)
        fig.canvas.draw_idle()

    # Auto-optimize checkbox
    ax_check = fig.add_axes([slider_x2, 0.28 - 2 * 0.04, 0.2, 0.03])
    check_auto = CheckButtons(ax_check, ['Auto-optimize gradient'], [params['auto_opt']])
    check_auto.labels[0].set_fontsize(9)
    set_grad_active(not params['auto_opt'])

    # Unit toggle button
    ax_unit_btn = fig.add_axes([slider_x2 + 0.22, 0.28 - 2 * 0.04, 0.06, 0.03])
    btn_unit = Button(ax_unit_btn, 'mm ↔ in', hovercolor='#ddd')
    btn_unit.label.set_fontsize(8)

    # Info text area (bottom right)
    ax_info = fig.add_axes([0.58, 0.02, 0.40, 0.12])
    ax_info.set_axis_off()
    info_text = ax_info.text(0, 1, '', fontsize=8.5, fontfamily='monospace',
                              verticalalignment='top', transform=ax_info.transAxes)

    # Beam info (bottom left)
    ax_beam = fig.add_axes([0.08, 0.02, 0.45, 0.06])
    ax_beam.set_axis_off()
    beam_text = ax_beam.text(0, 1, '', fontsize=8.5, fontfamily='monospace',
                              verticalalignment='top', transform=ax_beam.transAxes,
                              color='#666')

    # ---- Update function ----
    def update(val=None):
        ax.clear()

        bore_in    = s_bore.val
        lmag_in    = s_lmag.val
        spacing_in = s_spacing.val
        drift_in   = s_drift.val
        target_mm  = s_spot.val
        manual_g   = s_grad.val
        auto_opt   = params['auto_opt']

        bore_m  = bore_in * IN_TO_M
        lmag_m  = lmag_in * IN_TO_M
        gap_in  = spacing_in - lmag_in
        gap_m   = gap_in * IN_TO_M
        drift_m = drift_in * IN_TO_M
        target_m = target_mm * MM_TO_M

        if gap_m < 0.001:
            ax.text(0.5, 0.5, 'Spacing must be > magnet length',
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=14, color='red')
            fig.canvas.draw_idle()
            return

        # Get gradient
        if auto_opt:
            g1 = optimize_gradient(bore_m, lmag_m, gap_m, drift_m,
                                    params['energy'], params['x0'], params['xp0'],
                                    target_m, balance_tol=0.30)
        else:
            g1 = manual_g

        if g1 < 0.01:
            ax.text(0.5, 0.5,
                    'No balanced solution found.\n'
                    'Try larger bore, shorter drift, or relax spot size.',
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=12, color='red')
            info_text.set_text('')
            fig.canvas.draw_idle()
            return

        env = track_envelope(bore_m, lmag_m, gap_m, drift_m, g1,
                              params['energy'], params['x0'], params['xp0'])

        z = env['z']
        x = env['x']
        y_env = env['y']

        # Convert to display units
        z_d = np.array([units.m_to_display(v) for v in z])
        x_d = np.array([units.m_to_display(v) for v in x])
        y_d = np.array([units.m_to_display(v) for v in y_env])
        bore_d = units.m_to_display(bore_m)
        target_d = units.m_to_display(target_m)
        du = units.label

        # Region boundaries in display units
        q1_end_d   = units.m_to_display(env['q1_end'])
        q2_start_d = units.m_to_display(env['q2_start'])
        q2_end_d   = units.m_to_display(env['q2_end'])
        total_d    = units.m_to_display(env['total_length'])

        # ---- Plot ----
        # Quad shading
        ax.axvspan(0, q1_end_d, alpha=0.08, color='#4444cc', label='_')
        ax.axvspan(q2_start_d, q2_end_d, alpha=0.08, color='#4444cc', label='_')
        ax.text(q1_end_d / 2, bore_d * 0.92, 'Q1', ha='center', fontsize=9, color='#888', style='italic')
        ax.text((q2_start_d + q2_end_d) / 2, bore_d * 0.92, 'Q2', ha='center', fontsize=9, color='#888', style='italic')

        # Bore limits
        ax.axhline(bore_d, color='red', ls='--', lw=0.8, alpha=0.35, label=f'Bore ±{bore_in:.1f} in')
        ax.axhline(-bore_d, color='red', ls='--', lw=0.8, alpha=0.35)

        # Target spot lines
        ax.axhline(target_d, color='gray', ls=':', lw=0.6, alpha=0.5, label=f'Target ±{target_mm:.1f} mm')
        ax.axhline(-target_d, color='gray', ls=':', lw=0.6, alpha=0.5)

        # Beam envelopes (± symmetric)
        ax.plot(z_d, x_d, color='#2266bb', lw=2, label='x-plane')
        ax.plot(z_d, -x_d, color='#2266bb', lw=2)
        ax.plot(z_d, y_d, color='#cc5522', lw=2, label='y-plane')
        ax.plot(z_d, -y_d, color='#cc5522', lw=2)

        # Zero axis
        ax.axhline(0, color='black', lw=0.5, alpha=0.3)

        # Crossover markers
        for zc in env['crossovers_x']:
            zc_d = units.m_to_display(zc)
            ax.plot(zc_d, 0, 'o', color='#2266bb', ms=5, zorder=5)
        for zc in env['crossovers_y']:
            zc_d = units.m_to_display(zc)
            ax.plot(zc_d, 0, 's', color='#cc5522', ms=5, zorder=5)

        ax.set_xlabel(f'z ({du})', fontsize=10)
        ax.set_ylabel(f'r ({du})', fontsize=10)
        ax.set_xlim(-total_d * 0.02, total_d * 1.04)
        max_r = max(bore_d, np.max(np.abs(x_d)), np.max(np.abs(y_d))) * 1.15
        ax.set_ylim(-max_r, max_r)
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        ax.tick_params(labelsize=8)
        ax.set_title(f'Beam envelope through FD doublet  |  g = {g1:.2f} T/m', fontsize=10, pad=8)
        ax.grid(True, alpha=0.15)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(which='minor', length=2.5, color='gray')
        ax.grid(True, which='minor', alpha=0.06)

        # Clipping warning
        clipping = env['max_env_x'] > bore_m * 0.95 or env['max_env_y'] > bore_m * 0.95
        if clipping:
            ax.text(0.5, 0.95, '⚠ BEAM CLIPPED BY BORE',
                    ha='center', va='top', transform=ax.transAxes,
                    fontsize=11, color='red', fontweight='bold',
                    bbox=dict(boxstyle='round', facecolor='#fee', alpha=0.9))

        # ---- Compute outputs ----
        pole_tip_B = g1 * bore_m                      # [T]
        ampere_turns = pole_tip_B * bore_m / MU0       # [A·turns]

        final_x_mm = env['final_x'] * 1000
        final_y_mm = env['final_y'] * 1000
        avg_mm = (final_x_mm + final_y_mm) / 2
        imbalance = abs(final_x_mm - final_y_mm) / avg_mm * 100 if avg_mm > 0 else 0

        # Crossover strings
        cx_str = ', '.join([f"{units.m_to_display(zc):.1f}" for zc in env['crossovers_x']])
        cy_str = ', '.join([f"{units.m_to_display(zc):.1f}" for zc in env['crossovers_y']])
        if not cx_str: cx_str = 'none'
        if not cy_str: cy_str = 'none'

        # Format spot sizes in display units
        fx_d = units.m_to_display(env['final_x'])
        fy_d = units.m_to_display(env['final_y'])
        avg_d = (fx_d + fy_d) / 2
        spot_dec = 4 if units.use_inches else 2

        info_lines = [
            f"{'RESULTS':─^44}",
            f"  Gradient:        {g1:>8.2f}  T/m",
            f"  Pole-tip field:  {pole_tip_B*1000:>8.1f}  mT",
            f"  Ampere-turns:    {ampere_turns:>8.0f}  A·t",
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

        gap_d = units.to_display(gap_in)
        beam_lines = [
            f"Beam: {params['energy']} MeV proton  |  "
            f"Bρ = {Brho:.4f} T·m  |  "
            f"r₀ = {params['x0']*1000:.1f} mm  |  "
            f"θ₀ = {params['xp0']*1000:.0f} mrad  |  "
            f"Gap = {gap_d:.1f} {du}",
        ]
        beam_text.set_text('\n'.join(beam_lines))

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

        # Update slider labels (just the text, vals stay in original units)
        s_bore.label.set_text(f'Bore radius ({new_label})')
        s_lmag.label.set_text(f'Mag length ({new_label})')
        s_spacing.label.set_text(f'Q1–Q2 spacing ({new_label})')
        s_drift.label.set_text(f'Drift to target ({new_label})')
        s_spot.label.set_text(f'Target spot ({new_label})')

        # Update slider ranges and values for geometry (in ↔ mm)
        if units.use_inches:
            # Switched TO inches — sliders are already stored in inches
            s_bore.valmin, s_bore.valmax, s_bore.valstep = 1.5, 6.0, 0.1
            s_lmag.valmin, s_lmag.valmax, s_lmag.valstep = 2.0, 6.0, 0.1
            s_spacing.valmin, s_spacing.valmax, s_spacing.valstep = 3.0, 12.0, 0.1
            s_drift.valmin, s_drift.valmax, s_drift.valstep = 12.0, 24.0, 0.5
            s_spot.valmin, s_spot.valmax, s_spot.valstep = 0.1, 10.0, 0.1
        else:
            s_bore.valmin, s_bore.valmax, s_bore.valstep = 38.0, 152.0, 1.0
            s_lmag.valmin, s_lmag.valmax, s_lmag.valstep = 51.0, 152.0, 1.0
            s_spacing.valmin, s_spacing.valmax, s_spacing.valstep = 76.0, 305.0, 1.0
            s_drift.valmin, s_drift.valmax, s_drift.valstep = 305.0, 610.0, 5.0
            s_spot.valmin, s_spot.valmax, s_spot.valstep = 0.1, 10.0, 0.1

        # Note: slider .val stays in its original meaning (inches for geo, mm for spot)
        # The plot & info just display in the selected unit
        update()

    btn_unit.on_clicked(on_unit_toggle)

    # ---- Initial draw ----
    set_grad_active(not params['auto_opt'])
    update()

    print("=" * 60)
    print("  Quadrupole Doublet Design Tool")
    print("  1 MeV Proton | FD configuration")
    print("  Drag sliders or type values. Close window to exit.")
    print("=" * 60)
    print(f"  Bρ = {Brho:.4f} T·m")
    print(f"  Balance constraint: x/y spot within 30% of each other")
    print()

    plt.show()


if __name__ == '__main__':
    main()