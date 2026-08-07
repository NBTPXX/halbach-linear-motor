import json
import sys

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu
from scipy.signal import savgol_filter

# Linear 2D magnetic-vector-potential approximation.
# y is the travel direction and z points from the stator toward the magnets.
MU0 = 4 * np.pi * 1e-7
BR = 1.45
MU_R_IRON = 4000.0  # Linearized 35W300 working-region approximation.
HALBACH_WIDTHS = {"halbach66": 6.0, "halbach8": 8.0, "halbach10": 10.0, "halbach11": 11.0}
PITCH = 16.0
TOOTH = 5.0
PERIOD = 48.0
GAP = 1.0
H = 0.5
MAGNET_SUBCELL_SAMPLES = 8
COENERGY_SMOOTH_WINDOW = 9
COENERGY_SMOOTH_ORDER = 3
COIL_TURNS = 266
ACTIVE_STACK_M = 16e-3
IQ_PEAK_A = 4.0
TOOTH_CENTERS = np.array([2.5, 18.5, 34.5])
SINGLE_CORE = "--single-core" in sys.argv
CORE_OTHER_SIDE = "--core-other-side" in sys.argv
NO_SIDE = "--no-side" in sys.argv
THRUST_CURVE = "--thrust-curve" in sys.argv
CURVE_ONLY = "--curve-only" in sys.argv or THRUST_CURVE
ANIMATE = "--animate" in sys.argv or CURVE_ONLY
ANIMATION_FRAMES = int(sys.argv[sys.argv.index("--frames") + 1]) if "--frames" in sys.argv else 200
if ANIMATION_FRAMES < 4 or ANIMATION_FRAMES % 2:
    raise SystemExit("--frames must be an even integer of at least 4")
MAGNET_GROUPS = int(sys.argv[sys.argv.index("--magnet-groups") + 1]) if "--magnet-groups" in sys.argv else 5
CORE_CENTER_Y = 18.5
MAGNET_START = CORE_CENTER_Y - 12.0 * MAGNET_GROUPS
ANIMATION_TRAVEL = 24.0 if ANIMATE else 0.0
Y_MIN, Y_MAX = (CORE_CENTER_Y - 12.0 * MAGNET_GROUPS - 25.0 - ANIMATION_TRAVEL, CORE_CENTER_Y + 12.0 * MAGNET_GROUPS + 25.0 + ANIMATION_TRAVEL) if SINGLE_CORE else (0.0, PERIOD)
Z_MIN, Z_MAX = (-15.0, 30.0) if SINGLE_CORE and CORE_OTHER_SIDE else ((-25.0, 15.0) if SINGLE_CORE else (-20.0, 10.0))
Y_CELLS = int((Y_MAX - Y_MIN) / H)
Z_CELLS = int((Z_MAX - Z_MIN) / H)
Y = Y_MIN + (np.arange(Y_CELLS) + 0.5) * H
Z = Z_MIN + (np.arange(Z_CELLS) + 0.5) * H

# Closed DXF outline transformed from drawing coordinates into the solver frame:
# z_solver = y_dxf - 13 mm, placing the tooth tips at z = -1 mm.
CORE_DXF = np.array([
    (42.5, -5.0), (44.5, -5.0), (46.5, -7.0), (46.5, -19.5),
    (35.0, -19.5), (35.0, -17.5), (37.0, -17.5), (37.0, -15.5),
    (27.5, -15.5), (27.5, -17.5), (29.5, -17.5), (29.5, -19.5),
    (7.5, -19.5), (7.5, -17.5), (9.5, -17.5), (9.5, -15.5),
    (0.0, -15.5), (0.0, -17.5), (2.0, -17.5), (2.0, -19.5),
    (-9.5, -19.5), (-9.5, -7.0), (-7.5, -5.0), (-5.5, -5.0),
    (-5.5, -13.0), (0.0, -13.0), (0.0, -1.0), (5.0, -1.0),
    (5.0, -13.0), (16.0, -13.0), (16.0, -1.0), (21.0, -1.0),
    (21.0, -13.0), (32.0, -13.0), (32.0, -1.0), (37.0, -1.0),
    (37.0, -13.0), (42.5, -13.0),
])

# A 180 degree in-plane rotation puts the tooth tips at z = 4 mm: 1 mm above
# the 3 mm magnet thickness, while preserving the supplied asymmetric outline.
if CORE_OTHER_SIDE:
    CORE_DXF[:, 0] = 2 * CORE_CENTER_Y - CORE_DXF[:, 0]
    CORE_DXF[:, 1] = 3.0 - CORE_DXF[:, 1]


def inside_core(travel, height):
    inside = np.zeros_like(travel, dtype=bool)
    for start, end in zip(CORE_DXF, np.roll(CORE_DXF, -1, axis=0)):
        y1, z1 = start
        y2, z2 = end
        crosses = (z1 > height) != (z2 > height)
        crossing_y = (y2 - y1) * (height - z1) / (z2 - z1 + 1e-30) + y1
        inside ^= crosses & (travel < crossing_y)
    return inside


def material():
    mu_r = np.ones((Z_CELLS, Y_CELLS))
    travel, height = np.meshgrid(Y, Z)
    # Periodic force studies need neighboring cells; field renders use only the supplied outline.
    core = inside_core(travel, height)
    if not SINGLE_CORE:
        core |= inside_core(travel - PERIOD, height) | inside_core(travel + PERIOD, height)
    mu_r[core] = MU_R_IRON
    return mu_r


MU_R = material()
NU = 1.0 / (MU0 * MU_R)
IRON_AREA_MM2 = float(np.count_nonzero(MU_R > 1.0) * H * H)


def idx(j, i):
    return j * Y_CELLS + i


def build_operator():
    n = Y_CELLS * Z_CELLS
    matrix = lil_matrix((n, n))
    inv_h2 = 1.0 / (H * 1e-3) ** 2
    for j in range(Z_CELLS):
        for i in range(Y_CELLS):
            row = idx(j, i)
            # Open outer boundaries approximate a distant magnetic boundary.
            if j == 0 or j == Z_CELLS - 1 or (SINGLE_CORE and (i == 0 or i == Y_CELLS - 1)):
                matrix[row, row] = 1.0
                continue
            east_i = (i + 1) % Y_CELLS
            west_i = (i - 1) % Y_CELLS
            east = (NU[j, i] + NU[j, east_i]) * 0.5 * inv_h2
            west = (NU[j, i] + NU[j, west_i]) * 0.5 * inv_h2
            north = (NU[j, i] + NU[j + 1, i]) * 0.5 * inv_h2
            south = (NU[j, i] + NU[j - 1, i]) * 0.5 * inv_h2
            matrix[row, row] = -(east + west + north + south)
            matrix[row, idx(j, east_i)] = east
            matrix[row, idx(j, west_i)] = west
            matrix[row, idx(j + 1, i)] = north
            matrix[row, idx(j - 1, i)] = south
    return splu(matrix.tocsc())


SOLVER = build_operator()


def derivative_y(values):
    spacing = 2 * H * 1e-3
    if not SINGLE_CORE:
        return (np.roll(values, -1, axis=1) - np.roll(values, 1, axis=1)) / spacing
    derivative = np.zeros_like(values)
    derivative[:, 1:-1] = (values[:, 2:] - values[:, :-2]) / spacing
    derivative[:, 0] = (values[:, 1] - values[:, 0]) / (H * 1e-3)
    derivative[:, -1] = (values[:, -1] - values[:, -2]) / (H * 1e-3)
    return derivative


def magnet_br(position, skew_offset, array):
    br_y = np.zeros((Z_CELLS, Y_CELLS))
    br_z = np.zeros((Z_CELLS, Y_CELLS))
    in_thickness = (Z >= 0.0) & (Z < 3.0)
    # Average magnetization over each finite-volume cell. A binary assignment at
    # the cell centre makes a moving magnet edge jump by one full 0.5 mm column.
    subcell_y = Y[:, None] + H * (np.arange(MAGNET_SUBCELL_SAMPLES) + 0.5) / MAGNET_SUBCELL_SAMPLES - H / 2
    phase = subcell_y - position - skew_offset
    if not SINGLE_CORE:
        phase = np.mod(phase, 24.0)
        active_unit = np.ones_like(phase, dtype=bool)
    else:
        phase = phase - MAGNET_START
        active_unit = (phase >= 0.0) & (phase < 24.0 * MAGNET_GROUPS)
        phase = np.mod(phase, 24.0)
    if array in HALBACH_WIDTHS and not NO_SIDE:
        # One 24 mm pole pair: +Z, -Y, -Z, +Y. The normal/tangential widths sum to 12 mm.
        main_width = HALBACH_WIDTHS[array]
        br_z_line = np.where(active_unit & (phase < main_width), BR, np.where(active_unit & (phase >= 12.0) & (phase < 12.0 + main_width), -BR, 0.0)).mean(axis=1)
        br_y_line = np.where(active_unit & (phase >= main_width) & (phase < 12.0), -BR, np.where(active_unit & (phase >= 12.0 + main_width), BR, 0.0)).mean(axis=1)
    else:
        main_width = HALBACH_WIDTHS.get(array, 10.0)
        br_z_line = np.where(active_unit & (phase < main_width), BR, np.where(active_unit & (phase >= 12.0) & (phase < 12.0 + main_width), -BR, 0.0)).mean(axis=1)
        br_y_line = np.zeros_like(br_z_line)
    br_z[in_thickness, :] = br_z_line
    br_y[in_thickness, :] = br_y_line
    return br_y, br_z


def solve_fields(position, skew_offset, array):
    br_y, br_z = magnet_br(position, skew_offset, array)
    # curl(nu * Br) x-component for A_x formulation.
    rhs = derivative_y(NU * br_z)
    rhs -= (np.roll(NU * br_y, -1, axis=0) - np.roll(NU * br_y, 1, axis=0)) / (2 * H * 1e-3)
    rhs[0, :] = 0.0
    rhs[-1, :] = 0.0
    a = SOLVER.solve(rhs.ravel()).reshape((Z_CELLS, Y_CELLS))
    by = np.zeros_like(a)
    bz = np.zeros_like(a)
    by[1:-1, :] = (a[2:, :] - a[:-2, :]) / (2 * H * 1e-3)
    bz[:, :] = -derivative_y(a)
    return by, bz


def tooth_flux_linkages(bz):
    """Return flux linkages for the three 5 mm tooth faces in the air gap."""
    gap_row = np.argmin(np.abs(Z - 3.25))
    linkages = []
    for center in TOOTH_CENTERS:
        tooth = (Y >= center - TOOTH / 2) & (Y <= center + TOOTH / 2)
        flux_wb = np.trapezoid(bz[gap_row, tooth], Y[tooth] * 1e-3) * ACTIVE_STACK_M
        linkages.append(COIL_TURNS * flux_wb)
    return np.asarray(linkages)


def slice_energy(position, skew_offset, array):
    by, bz = solve_fields(position, skew_offset, array)
    # Coenergy density in the linearized material model, per metre of x width.
    density = (by * by + bz * bz) / (2 * MU0 * MU_R)
    return density.sum() * (H * 1e-3) ** 2


positions = np.arange(0.0, PERIOD, 0.5)
array = "10to2"
if "--array" in sys.argv:
    array = sys.argv[sys.argv.index("--array") + 1]
# 26.6 degree skew creates an 8 mm shift across the 16 mm active stack.
angle_deg = 0.0
if "--angle" in sys.argv:
    angle_deg = float(sys.argv[sys.argv.index("--angle") + 1])
skew_shift = 16.0 * np.tan(np.deg2rad(angle_deg))
# Resolve the 4 mm dominant cogging period with at most 0.5 mm spacing.
# Five fixed slices alias when the total skew is an integer multiple of 4 mm.
skew_slices = max(5, int(np.ceil(abs(skew_shift) / 0.5)) + 1)
skew_offsets = np.linspace(-skew_shift / 2, skew_shift / 2, skew_slices)

if "--render" in sys.argv:
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle

    output_path = sys.argv[sys.argv.index("--render") + 1]
    render_offset = skew_shift / 2 if SINGLE_CORE and abs(skew_shift) > 1e-9 else 0.0
    by, bz = solve_fields(0.0, render_offset, array)
    magnitude_mt = np.hypot(by, bz) * 1e3
    if SINGLE_CORE:
        view_y, view_by, view_bz, view_magnitude_mt = Y, by, bz, magnitude_mt
    else:
        view_y = np.r_[Y - PERIOD, Y, Y + PERIOD]
        view_by = np.tile(by, (1, 3))
        view_bz = np.tile(bz, (1, 3))
        view_magnitude_mt = np.tile(magnitude_mt, (1, 3))
    figure, axis = plt.subplots(figsize=(14, 7), constrained_layout=True)
    image = axis.imshow(np.log10(np.clip(view_magnitude_mt, 1, 2000)), origin="lower", extent=(view_y[0] - H / 2, view_y[-1] + H / 2, Z_MIN, Z_MAX), aspect="equal", cmap="turbo")
    axis.streamplot(view_y, Z, view_by, view_bz, density=2.4, color="white", linewidth=0.55, arrowsize=0.6)
    profile = CORE_DXF
    axis.plot(np.r_[profile[:, 0], profile[0, 0]], np.r_[profile[:, 1], profile[0, 1]], color="white", linewidth=1.35, zorder=5)
    if array in HALBACH_WIDTHS and not NO_SIDE:
        main_width = HALBACH_WIDTHS[array]
        side_width = 12.0 - main_width
        blocks = [(0, main_width, "+Z"), (main_width, side_width, "-Y"), (12, main_width, "-Z"), (12 + main_width, side_width, "+Y")]
    else:
        main_width = HALBACH_WIDTHS.get(array, 10.0)
        blocks = [(0, main_width, "+Z"), (12, main_width, "-Z")]
    group_count = MAGNET_GROUPS if SINGLE_CORE else 1
    for group in range(group_count):
        group_start = (MAGNET_START if SINGLE_CORE else 0.0) + group * 24.0 + render_offset
        for start, width, label in blocks:
            axis.add_patch(Rectangle((group_start + start, 0), width, 3, facecolor="none", edgecolor="white", linewidth=1.35, zorder=6))
            axis.text(group_start + start + width / 2, 1.5, label, ha="center", va="center", color="white", fontsize=8, weight="bold", zorder=7)
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label("log10 |B| (mT)")
    model_name = "single DXF iron core with finite magnet array" if SINGLE_CORE else "periodic DXF iron-core cells"
    if CORE_OTHER_SIDE:
        model_name += ", core rotated 180 degrees above magnets"
    if NO_SIDE:
        model_name += ", no side magnets"
    if abs(angle_deg) > 1e-9:
        model_name += f", {angle_deg:g} degree skew edge slice"
    axis.set(title=f"N52 {array} Halbach | {model_name}", xlabel="Travel Y (mm)", ylabel="Height Z (mm)", xlim=(Y_MIN, Y_MAX), ylim=(Z_MIN, Z_MAX))
    figure.savefig(output_path, dpi=180)
    print(json.dumps({"field_map": output_path, "array": array, "iron_area_mm2": round(IRON_AREA_MM2, 3)}, ensure_ascii=False))
    sys.exit(0)

if ANIMATE:
    from matplotlib import animation
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle

    output_path = sys.argv[sys.argv.index("--animate") + 1] if "--animate" in sys.argv else None
    sample_positions = np.linspace(-24.0, 24.0, ANIMATION_FRAMES // 2)
    core_positions = np.concatenate((sample_positions, sample_positions[::-1]))
    magnitudes_mt = []
    energies = []
    phase_flux_linkages = []
    visual_skew_offset = skew_shift / 2 if abs(skew_shift) > 1e-9 else 0.0
    force_skew_offsets = np.linspace(-skew_shift / 2, skew_shift / 2, 3) if abs(skew_shift) > 1e-9 else np.array([0.0])
    for core_position in sample_positions:
        # Translating the magnet source by -core_position is equivalent to
        # translating the core by core_position in the fixed-magnet view.
        by, bz = solve_fields(-core_position, visual_skew_offset, array)
        magnitudes_mt.append(np.hypot(by, bz) * 1e3)
        slice_energies = []
        slice_linkages = []
        for offset in force_skew_offsets:
            if abs(offset - visual_skew_offset) < 1e-9:
                slice_by, slice_bz = by, bz
            else:
                slice_by, slice_bz = solve_fields(-core_position, offset, array)
            density = (slice_by * slice_by + slice_bz * slice_bz) / (2 * MU0 * MU_R)
            slice_energies.append(density.sum() * (H * 1e-3) ** 2)
            slice_linkages.append(tooth_flux_linkages(slice_bz))
        energies.append(float(np.mean(slice_energies)) * 16e-3)
        phase_flux_linkages.append(np.mean(slice_linkages, axis=0))
    sample_step_m = (sample_positions[1] - sample_positions[0]) * 1e-3
    # Differentiate a local polynomial fit to coenergy rather than raw cell-level
    # energy samples. This preserves the resolved 4 mm cogging content while
    # rejecting residual sub-cell boundary noise before it becomes a force spike.
    force_samples = -savgol_filter(
        np.asarray(energies),
        COENERGY_SMOOTH_WINDOW,
        COENERGY_SMOOTH_ORDER,
        deriv=1,
        delta=sample_step_m,
        mode="interp",
    )
    phase_flux_linkages = np.asarray(phase_flux_linkages)
    electrical_angle = 2 * np.pi * sample_positions[:, None] / 24.0
    # Identify each phase d-axis from its computed fundamental flux linkage.
    # This keeps the q-axis current orthogonal to the actual finite-array field.
    d_axis = np.angle(np.sum(phase_flux_linkages * np.exp(-1j * electrical_angle), axis=0))
    current_angle = electrical_angle + d_axis
    # SVPWM adds common-mode voltage, while delta branch currents remain balanced
    # fundamental q-axis currents. Iq here is the winding-branch peak current.
    phase_currents = -IQ_PEAK_A * np.sin(current_angle)
    linkage_gradient = np.gradient(phase_flux_linkages, sample_positions * 1e-3, axis=0)
    thrust_samples = np.sum(phase_currents * linkage_gradient, axis=1)
    forces = np.concatenate((force_samples, force_samples[::-1]))
    magnitudes_mt = magnitudes_mt + list(reversed(magnitudes_mt))
    max_force = max(float(np.max(np.abs(forces))), 1e-9)

    if array in HALBACH_WIDTHS and not NO_SIDE:
        main_width = HALBACH_WIDTHS[array]
        side_width = 12.0 - main_width
        blocks = [(0, main_width, "+Z"), (main_width, side_width, "-Y"), (12, main_width, "-Z"), (12 + main_width, side_width, "+Y")]
    else:
        main_width = HALBACH_WIDTHS.get(array, 10.0)
        blocks = [(0, main_width, "+Z"), (12, main_width, "-Z")]

    figure, axis = plt.subplots(figsize=(13, 6), constrained_layout=True)
    display_min = MAGNET_START - 25.0
    display_max = MAGNET_START + MAGNET_GROUPS * 24.0 + 25.0

    def draw_frame(frame_index):
        core_position = core_positions[frame_index]
        force = forces[frame_index]
        axis.clear()
        magnitude_mt = magnitudes_mt[frame_index]
        axis.imshow(
            np.log10(np.clip(magnitude_mt, 1, 2000)),
            origin="lower",
            extent=(Y[0] + core_position - H / 2, Y[-1] + core_position + H / 2, Z_MIN, Z_MAX),
            aspect="equal",
            cmap="turbo",
            vmin=0,
            vmax=np.log10(2000),
        )
        profile = CORE_DXF + np.array([core_position, 0.0])
        axis.plot(np.r_[profile[:, 0], profile[0, 0]], np.r_[profile[:, 1], profile[0, 1]], color="white", linewidth=1.5, zorder=5)
        for group in range(MAGNET_GROUPS):
            group_start = MAGNET_START + group * 24.0
            for start, width, label in blocks:
                axis.add_patch(Rectangle((group_start + start, 0), width, 3, facecolor="none", edgecolor="white", linewidth=1.25, zorder=6))
                axis.text(group_start + start + width / 2, 1.5, label, ha="center", va="center", color="white", fontsize=7, weight="bold", zorder=7)
        arrow_length = 5.0 + 18.0 * abs(force) / max_force
        core_center = CORE_CENTER_Y + core_position
        arrow_start = core_center - arrow_length / 2 if force >= 0 else core_center + arrow_length / 2
        axis.annotate("", xy=(arrow_start + np.sign(force if force else 1.0) * arrow_length, 13.0), xytext=(arrow_start, 13.0), arrowprops={"arrowstyle": "->", "color": "white", "lw": 2.8}, zorder=8)
        axis.text(core_center, 15.0, f"Position {core_position:+.0f} mm\nCogging force {force:+.3f} N", ha="center", va="bottom", color="white", fontsize=11, weight="bold", zorder=8, bbox={"boxstyle": "round,pad=0.35", "facecolor": "#111827", "edgecolor": "white", "alpha": 0.9})
        axis.set(
            title="N52 6:6 Halbach | Core motion and cogging force",
            xlabel="Travel Y (mm)",
            ylabel="Height Z (mm)",
            xlim=(display_min, display_max),
            ylim=(Z_MIN, Z_MAX),
        )

    if "--animation-preview" in sys.argv:
        preview_path = sys.argv[sys.argv.index("--animation-preview") + 1]
        draw_frame(0)
        figure.savefig(preview_path, dpi=130)
        print(json.dumps({"animation_preview": preview_path, "frame": 0}, ensure_ascii=False))
        sys.exit(0)

    curve_path = None
    thrust_path = None
    if "--curve" in sys.argv:
        curve_path = sys.argv[sys.argv.index("--curve") + 1]
        curve_figure, curve_axis = plt.subplots(figsize=(11, 4.5), constrained_layout=True)
        curve_axis.plot(sample_positions, force_samples, color="#2563eb", linewidth=2.2, label="Total cogging force")
        curve_axis.fill_between(sample_positions, force_samples, 0.0, color="#2563eb", alpha=0.12)
        force_fft = np.fft.rfft(force_samples - force_samples.mean())
        component_indexes = np.argsort(np.abs(force_fft[1:]))[-3:] + 1
        component_indexes = component_indexes[np.argsort(np.abs(force_fft[component_indexes]))[::-1]]
        display_positions = np.linspace(sample_positions[0], sample_positions[-1], 1000)
        display_indexes = (display_positions - sample_positions[0]) / (sample_positions[1] - sample_positions[0])
        sample_step = sample_positions[1] - sample_positions[0]
        component_colors = ("#dc2626", "#d97706", "#7c3aed")
        for component_index, color in zip(component_indexes, component_colors):
            amplitude = 2 * abs(force_fft[component_index]) / len(force_samples)
            phase = np.angle(force_fft[component_index])
            component = amplitude * np.cos(2 * np.pi * component_index * display_indexes / len(force_samples) + phase)
            period_mm = sample_step * len(force_samples) / component_index
            curve_axis.plot(display_positions, component, "--", color=color, linewidth=1.35, label=f"FFT {period_mm:.2f} mm, A={amplitude:.3f} N")
        curve_axis.axhline(0.0, color="#64748b", linewidth=1)
        curve_axis.set(title="Cogging force with leading FFT components", xlabel="Core position x (mm)", ylabel="Cogging force F (N)", xlim=(-24, 24))
        curve_axis.grid(alpha=0.25)
        curve_axis.legend(loc="upper right", fontsize=8)
        curve_figure.savefig(curve_path, dpi=160)
        plt.close(curve_figure)

    if THRUST_CURVE:
        thrust_path = sys.argv[sys.argv.index("--thrust-curve") + 1]
        thrust_figure, thrust_axis = plt.subplots(figsize=(11, 4.5), constrained_layout=True)
        mean_thrust = float(np.mean(thrust_samples))
        thrust_axis.plot(sample_positions, thrust_samples, color="#059669", linewidth=2.2, label="FOC thrust")
        thrust_axis.fill_between(sample_positions, thrust_samples, mean_thrust, color="#059669", alpha=0.12)
        thrust_axis.axhline(mean_thrust, color="#475569", linestyle="--", linewidth=1.2, label=f"Mean = {mean_thrust:.3f} N")
        thrust_fft = np.fft.rfft(thrust_samples - mean_thrust)
        component_indexes = np.argsort(np.abs(thrust_fft[1:]))[-3:] + 1
        component_indexes = component_indexes[np.argsort(np.abs(thrust_fft[component_indexes]))[::-1]]
        display_positions = np.linspace(sample_positions[0], sample_positions[-1], 1000)
        display_indexes = (display_positions - sample_positions[0]) / (sample_positions[1] - sample_positions[0])
        sample_step = sample_positions[1] - sample_positions[0]
        component_colors = ("#dc2626", "#d97706", "#7c3aed")
        for component_index, color in zip(component_indexes, component_colors):
            amplitude = 2 * abs(thrust_fft[component_index]) / len(thrust_samples)
            phase = np.angle(thrust_fft[component_index])
            component = mean_thrust + amplitude * np.cos(2 * np.pi * component_index * display_indexes / len(thrust_samples) + phase)
            period_mm = sample_step * len(thrust_samples) / component_index
            thrust_axis.plot(display_positions, component, "--", color=color, linewidth=1.35, label=f"FFT {period_mm:.2f} mm, A={amplitude:.3f} N")
        thrust_axis.set(
            title="SVPWM delta FOC thrust with leading FFT components at Iq = 4 A",
            xlabel="Core position x (mm)",
            ylabel="Thrust F (N)",
            xlim=(-24, 24),
        )
        thrust_axis.grid(alpha=0.25)
        thrust_axis.legend(loc="upper right", fontsize=8)
        thrust_figure.savefig(thrust_path, dpi=160)
        plt.close(thrust_figure)

    if CURVE_ONLY:
        print(json.dumps({
            "curve": curve_path,
            "thrust_curve": thrust_path,
            "samples": len(sample_positions),
            "force_peak_to_peak_N": round(float(forces.max() - forces.min()), 4),
            "thrust_mean_N": round(float(np.mean(thrust_samples)), 4),
            "thrust_peak_to_peak_N": round(float(np.ptp(thrust_samples)), 4),
            "iq_peak_A": IQ_PEAK_A,
        }, ensure_ascii=False))
        sys.exit(0)

    movie = animation.FuncAnimation(figure, draw_frame, frames=len(core_positions), interval=120)
    movie.save(output_path, writer=animation.PillowWriter(fps=15), dpi=110)
    print(json.dumps({"animation": output_path, "curve": curve_path if "--curve" in sys.argv else None, "frames": len(core_positions), "travel_mm": [-24, 24, -24], "force_peak_to_peak_N": round(float(forces.max() - forces.min()), 4)}, ensure_ascii=False))
    sys.exit(0)

energies = []
for position in positions:
    slice_energies = np.array([slice_energy(position, offset, array) for offset in skew_offsets])
    if abs(skew_shift) > 1e-9:
        per_width = np.trapezoid(slice_energies, skew_offsets) / abs(skew_shift)
    else:
        per_width = slice_energies[0]
    energies.append(per_width * 16e-3)
energies = np.array(energies)
forces = -(np.roll(energies, -1) - np.roll(energies, 1)) / (1e-3)

harmonics = []
for order in range(1, 25):
    phase = 2 * np.pi * order * np.arange(len(forces)) / len(forces)
    amplitude = 2 * np.hypot(np.sum(forces * np.cos(phase)), np.sum(forces * np.sin(phase))) / len(forces)
    if amplitude > 0.001:
        harmonics.append({"period_mm": round(PERIOD / order, 3), "amplitude_N": round(float(amplitude), 4)})

print(json.dumps({
    "mean_coenergy_mJ": round(float(energies.mean() * 1e3), 4),
    "cogging_peak_to_peak_N": round(float(forces.max() - forces.min()), 4),
    "cogging_peak_N": round(float(np.max(np.abs(forces))), 4),
    "cogging_rms_N": round(float(np.sqrt(np.mean(forces ** 2))), 4),
    "skew_angle_deg": angle_deg,
    "skew_shift_mm": round(float(skew_shift), 4),
    "skew_slices": skew_slices,
    "array": array,
    "iron_area_mm2": round(IRON_AREA_MM2, 3),
    "harmonics": harmonics,
}, ensure_ascii=False, indent=2))
