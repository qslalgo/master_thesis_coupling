import os
import glob
import argparse
import numpy as np
import openmc

from model_builder import build_model
from power_extraction import build_qprime_pin_for_solver_by_assembly
from dual_velocity_solver import Dual_Velocity_Assembly_Solver
from cases import cases
from model_reconstruct import (
    reconstruct_model,
    build_fuel_slices,
    build_clad_slices,
    build_coolant_slices,
    create_axial_planes,
)
from iapws import IAPWS97


# ============================================================
# Utilities
# ============================================================

def water_density_iapws(T_K, P_MPa=0.1013):
    T = np.asarray(T_K, dtype=float)
    rho = np.zeros_like(T)
    for i, Ti in np.ndenumerate(T):
        rho[i] = IAPWS97(T=float(Ti), P=P_MPa).rho   # kg/m³
    return rho / 1000.0                                # g/cm³


def flatten_cell_dict(cell_dict):
    return [c for aid in sorted(cell_dict) for c in cell_dict[aid]]


def get_mesh_path(mesh_dir, config, assembly_type):
    return os.path.join(mesh_dir, config.mesh_map[assembly_type])


def is_th_result_degenerate(res, config):
    return np.allclose(res["T_inner_z_K"], 273.15, atol=1e-3)


# ============================================================
# Logging
# ============================================================

def log_iteration(log_file, it, th_results):
    log_file.write(f"\n=== Iteration {it} ===\n")
    for assy_id in sorted(th_results):
        r = th_results[assy_id]
        log_file.write(
            f"Assembly {assy_id}: "
            f"T_in={r['T_inner_z_K'][0]:.2f}, "
            f"T_out={r['T_inner_z_K'][-1]:.2f}, "
            f"rho_in={r['rho_inner_z'][0]:.5f}, "
            f"rho_out={r['rho_inner_z'][-1]:.5f}, "
            f"P_gen={r['P_gen']:.1f}, "
            f"P_in={r['P_to_inner']:.1f}, "
            f"P_out={r['P_to_outer']:.1f}, "
            f"P_mcp={r['P_mcpDT_in']:.1f}, "
            f"mdot_in={r['mdot_in']:.6e}, "
            f"mdot_out={r['mdot_out']:.6e}, "
            f"exact_f_Re_inner={r['exact_f_Re_inner']:.6e}, "
            f"A_flow_inner={r['A_flow_inner']:.6e}, "
            f"A_flow_outer={r['A_flow_outer']:.6e}, "
            f"D_h_inner={r['D_h_inner']:.6e}, "
            f"D_h_outer={r['D_h_outer']:.6e}, "
            f"A_holes_bottom={r['A_holes_bottom']:.6e}, "
            f"K_bottom={r['K_bottom']:.6e}, "
            f"K_top={r['K_top']:.6e}\n"
        )
    log_file.flush()


# ============================================================
# OpenMC runner
# ============================================================

def run_openmc_in_dir(model, workdir, PlotColors=None):
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    model.export_to_xml()
    print(model.geometry.bounding_box)
    openmc.run()
    sp_path = sorted(glob.glob("statepoint.*.h5"))[-1]
    return os.path.join(workdir, sp_path)


# ============================================================
# TH runner
# ============================================================

def run_th_case(config, mesh_path, qprime_data, case_name, outdir):
    solver = Dual_Velocity_Assembly_Solver(config, mesh_path)
    solver.assign_qprime_by_centers(qprime_data["values"], qprime_data["centers"])
    solver.run_old(case_name=case_name, out_dir=outdir)

    T_inner_z_K = np.asarray(solver.plot_T_bulk_in)  + 273.15
    T_outer_z_K = np.asarray(solver.plot_T_bulk_out) + 273.15
    T_fuel_z_K  = np.asarray(solver.plot_T_fuel)     + 273.15
    T_clad_z_K  = 0.5 * (T_fuel_z_K + T_inner_z_K)

    rho_inner_z = water_density_iapws(T_inner_z_K)
    rho_outer_z = water_density_iapws(T_outer_z_K)

    P_gen      = np.sum(solver.qprime_pin) * solver.dz
    P_to_inner = np.sum(solver.plot_heat_in)  * solver.dz
    P_to_outer = np.sum(solver.plot_heat_out) * solver.dz

    dT         = np.diff(solver.plot_T_bulk_in, prepend=solver.plot_T_bulk_in[0])
    P_mcpDT_in = np.sum(solver.plot_mdot_in * solver.plot_cp_in * dT)

    print(f"  P_gen={P_gen:.1f} W  P_inner={P_to_inner:.1f} W  "
          f"P_outer={P_to_outer:.1f} W  P_mcp={P_mcpDT_in:.1f} W")

    return {
        "qprime_pin":       solver.qprime_pin,
        "T_inner_z_K":      T_inner_z_K,
        "T_outer_z_K":      T_outer_z_K,
        "T_fuel_z_K":       T_fuel_z_K,
        "T_clad_z_K":       T_clad_z_K,
        "rho_inner_z":      rho_inner_z,
        "rho_outer_z":      rho_outer_z,
        "P_gen":            P_gen,
        "P_to_inner":       P_to_inner,
        "P_to_outer":       P_to_outer,
        "P_mcpDT_in":       P_mcpDT_in,
        "mdot_in":          solver.plot_mdot_in,
        "mdot_out":         solver.plot_mdot_out,
        "exact_f_Re_inner": solver.exact_f_Re_inner,
        "A_flow_inner":     solver.A_flow_inner,
        "A_flow_outer":     solver.A_flow_outer,
        "D_h_inner":        solver.D_h_inner,
        "D_h_outer":        solver.D_h_outer,
        "A_holes_bottom":   solver.A_holes_bottom,
        "K_bottom":         solver.K_bottom,
        "K_top":            solver.K_top,
    }


# ============================================================
# Particles schedule
# ============================================================

_PARTICLES = [
    (0,  0,    1_000,  100,  20),
    (1,  3,    3_000,  300,  80),
    (4,  6,    6_000,  600, 100),
    (7,  10,   8_000,  800, 100),
    (11, 9999, 10_000, 1_000, 200),
]

def particles_for_iteration(it):
    for lo, hi, p, b, i in _PARTICLES:
        if lo <= it <= hi:
            return p, b, i
    return 10_000, 1_000, 200


# ============================================================
# Main feedback loop
# ============================================================

def full_feedback_loop(base_dir, mesh_path, case_name, config,
                       max_iter=10, alpha=1.0, log_file=None):
    os.chdir(base_dir)

    # ── build model ──────────────────────────────────────────
    model, geometry_data, outer_coolant_cells, PlotColors = build_model(config)
    universe_0 = model.geometry.root_universe

    z_planes = create_axial_planes(config.z_min, config.z_max, config.n_z)

    fuel_cells_by_assembly          = geometry_data["fuel_cells"]
    clad_cells_by_assembly          = geometry_data["clad_cells"]
    inner_coolant_cells_by_assembly = geometry_data["inner_coolant_cells"]
    assembly_types                  = geometry_data["assembly_types_by_id_openmc"]

    water_below_cells_by_assembly   = geometry_data["water_below_assembly_cells"]
    water_above_cells_by_assembly   = geometry_data["water_above_assembly_cells"]
    top_grid_cells_by_assembly      = geometry_data["top_grid_cells"]
    bottom_grid_cells_by_assembly   = geometry_data["bottom_grid_cells"]

    # ── axial slicing ─────────────────────────────────────────
    fuel_layers_by_assembly          = {}
    clad_layers_by_assembly          = {}
    inner_coolant_layers_by_assembly = {}
    inner_coolant_materials          = {}

    for aid, cells in fuel_cells_by_assembly.items():
        fuel_layers_by_assembly[aid] = build_fuel_slices(
            cells, universe_0, z_planes, config.n_z)

    for aid, cells in clad_cells_by_assembly.items():
        layers, lower, upper = build_clad_slices(
            cells, universe_0, z_planes, config.n_z)
        clad_layers_by_assembly[aid] = {"layers": layers, "lower": lower, "upper": upper}

    for aid, cells in inner_coolant_cells_by_assembly.items():
        layers, lower_cell, upper_cell = build_coolant_slices(
            cells, universe_0, z_planes, config.n_z)
        inner_coolant_layers_by_assembly[aid] = {
            "layers": layers, "lower": lower_cell, "upper": upper_cell}

        base_mat = cells[0].fill
        inner_coolant_materials[aid] = []
        for k in range(config.n_z):
            m = base_mat.clone(); m.name = f"coolant_a{aid}_z{k}"
            inner_coolant_materials[aid].append(m)

    for aid, data in inner_coolant_layers_by_assembly.items():
        for k, cells in enumerate(data["layers"]):
            for c in cells:
                c.fill = inner_coolant_materials[aid][k]

    # boundary materials
    inner_coolant_boundary_materials = {}
    for aid, cells in inner_coolant_cells_by_assembly.items():
        base_mat = cells[0].fill
        lo = base_mat.clone(); lo.name = f"coolant_a{aid}_lower"
        hi = base_mat.clone(); hi.name = f"coolant_a{aid}_upper"
        inner_coolant_boundary_materials[aid] = {"lower": lo, "upper": hi}

    for aid, data in inner_coolant_layers_by_assembly.items():
        data["lower"].fill = inner_coolant_boundary_materials[aid]["lower"]
        data["upper"].fill = inner_coolant_boundary_materials[aid]["upper"]

    model.materials.extend(m for mats in inner_coolant_materials.values() for m in mats)
    model.materials.extend(m for d in inner_coolant_boundary_materials.values()
                           for m in d.values())

    # ── outer coolant ─────────────────────────────────────────
    outer_layers, outer_lower, outer_upper = build_coolant_slices(
        outer_coolant_cells, universe_0, z_planes, config.n_z)
    outer_coolant_data = {"layers": outer_layers, "lower": outer_lower, "upper": outer_upper}

    base_outer = outer_coolant_cells[0].fill
    outer_mats = {
        "layers": {},
        "lower":  (lambda m: setattr(m, "name", "outer_coolant_lower") or m)(base_outer.clone()),
        "upper":  (lambda m: setattr(m, "name", "outer_coolant_upper") or m)(base_outer.clone()),
    }
    for k in range(config.n_z):
        m = base_outer.clone(); m.name = f"outer_coolant_z{k}"
        outer_mats["layers"][k] = m
        for c in outer_layers[k]:
            c.fill = m

    outer_lower.fill = outer_mats["lower"]
    outer_upper.fill = outer_mats["upper"]
    model.materials.extend(outer_mats["layers"].values())
    model.materials.append(outer_mats["lower"])
    model.materials.append(outer_mats["upper"])

    # ── water above / below ───────────────────────────────────
    water_above_mats = {}
    water_below_mats = {}
    for aid, cells in water_above_cells_by_assembly.items():
        water_above_mats[aid] = []
        for cell in cells:
            m = cell.fill.clone(); cell.fill = m
            water_above_mats[aid].append(m)
    for aid, cells in water_below_cells_by_assembly.items():
        water_below_mats[aid] = []
        for cell in cells:
            m = cell.fill.clone(); cell.fill = m
            water_below_mats[aid].append(m)

    model.materials.extend(m for mats in water_above_mats.values() for m in mats)
    model.materials.extend(m for mats in water_below_mats.values() for m in mats)

    # ── initial TH state ──────────────────────────────────────
    T_init = getattr(config, "T_inlet_K", 295.65)

    T_fuel    = {aid: np.full(config.n_z, T_init) for aid in fuel_cells_by_assembly}
    T_clad    = {aid: np.full(config.n_z, T_init) for aid in fuel_cells_by_assembly}
    T_inner   = {aid: np.full(config.n_z, T_init) for aid in fuel_cells_by_assembly}
    rho_inner = {aid: np.full(config.n_z, water_density_iapws(T_init))
                 for aid in fuel_cells_by_assembly}
    T_outer   = np.full(config.n_z, T_init)
    rho_outer = np.full(config.n_z, water_density_iapws(T_init))

    # ── convergence bookkeeping (outside loop!) ───────────────
    eps             = 1e-8
    tol_power_total = 1e-2
    tol_power_assy  = 2e-2
    P_total_prev    = config.P_core_W * 1.44
    P_per_assy_prev = np.full(len(assembly_types),
                              config.P_core_W * 1.44 / len(assembly_types))

    # ============================================================
    # Iteration loop
    # ============================================================
    for it in range(max_iter):
        print(f"\n{'='*16} ITERATION {it} {'='*16}")

        # ── apply temperatures / densities to OpenMC ──────────
        for aid, layers in fuel_layers_by_assembly.items():
            for k in range(config.n_z):
                for c in layers[k]:
                    c.temperature = float(T_fuel[aid][k])

        for aid, data in clad_layers_by_assembly.items():
            for k, cells in enumerate(data["layers"]):
                for c in cells:
                    c.temperature = float(T_clad[aid][k])
            data["lower"].temperature = float(T_clad[aid][0])
            data["upper"].temperature = float(T_clad[aid][-1])

        for aid, data in inner_coolant_layers_by_assembly.items():
            for k, cells in enumerate(data["layers"]):
                for c in cells:
                    c.temperature = float(T_inner[aid][k])
                inner_coolant_materials[aid][k].set_density(
                    "g/cm3", float(rho_inner[aid][k]))

            T_bot = float(T_inner[aid][0]);  rho_bot = float(rho_inner[aid][0])
            T_top = float(T_inner[aid][-1]); rho_top = float(rho_inner[aid][-1])

            data["lower"].temperature = T_bot
            inner_coolant_boundary_materials[aid]["lower"].set_density("g/cm3", rho_bot)
            data["upper"].temperature = T_top
            inner_coolant_boundary_materials[aid]["upper"].set_density("g/cm3", rho_top)

            for cell in water_above_cells_by_assembly[aid]:
                cell.temperature = T_top
                cell.fill.set_density("g/cm3", rho_top)
            for cell in water_below_cells_by_assembly[aid]:
                cell.temperature = T_bot
                cell.fill.set_density("g/cm3", rho_bot)
            for cell in top_grid_cells_by_assembly[aid]:
                cell.temperature = T_top
            for cell in bottom_grid_cells_by_assembly[aid]:
                cell.temperature = T_bot

        for k in range(config.n_z):
            for c in outer_layers[k]:
                c.temperature = float(T_outer[k])
            outer_mats["layers"][k].set_density("g/cm3", float(rho_outer[k]))

        outer_coolant_data["lower"].temperature = float(T_outer[0])
        outer_mats["lower"].set_density("g/cm3", float(rho_outer[0]))
        outer_coolant_data["upper"].temperature = float(T_outer[-1])
        outer_mats["upper"].set_density("g/cm3", float(rho_outer[-1]))

        model.materials = openmc.Materials(model.geometry.get_all_materials().values())

        # ── OpenMC run ────────────────────────────────────────
        particles, batches, inactive = particles_for_iteration(it)
        model.settings.particles = particles
        model.settings.batches   = batches
        model.settings.inactive  = inactive
        print(f"OpenMC: particles={particles}, batches={batches}, inactive={inactive}")

        omc_dir = os.path.join(base_dir, f"case_{case_name}", f"it_{it:02d}", "openmc")
        sp_path = run_openmc_in_dir(model, omc_dir, PlotColors=PlotColors)

        qprime_by_assembly = build_qprime_pin_for_solver_by_assembly(sp_path, config)

        # ── TH solve ─────────────────────────────────────────
        th_dir = os.path.join(base_dir, f"case_{case_name}", f"it_{it:02d}", "th")
        os.makedirs(th_dir, exist_ok=True)

        th_results = {}
        for aid, assy_type in assembly_types.items():
            print(f"  TH: assembly {aid} ({assy_type})")
            th_results[aid] = run_th_case(
                config,
                get_mesh_path(mesh_path, config, assy_type),
                qprime_by_assembly[aid],
                case_name=f"{case_name}_it{it:02d}_assy{aid}",
                outdir=os.path.join(th_dir, f"assy_{aid}"),
            )

        # ── relaxed update ────────────────────────────────────
        max_dT_fuel_rel  = 0.0
        max_dT_inner_rel = 0.0

        for aid, res in th_results.items():
            if is_th_result_degenerate(res, config):
                print(f"  ⚠️  Assembly {aid}: TH degenerate — reset to T_inlet_K")
                T_fuel[aid]    = np.full(config.n_z, config.T_inlet_K)
                T_clad[aid]    = np.full(config.n_z, config.T_inlet_K)
                T_inner[aid]   = np.full(config.n_z, config.T_inlet_K)
                rho_inner[aid] = np.full(config.n_z, water_density_iapws(config.T_inlet_K))
                continue

            T_fuel_new    = alpha * res["T_fuel_z_K"]  + (1 - alpha) * T_fuel[aid]
            T_clad_new    = alpha * res["T_clad_z_K"]  + (1 - alpha) * T_clad[aid]
            T_inner_new   = alpha * res["T_inner_z_K"] + (1 - alpha) * T_inner[aid]
            rho_inner_new = water_density_iapws(T_inner_new)

            max_dT_fuel_rel = max(max_dT_fuel_rel, np.max(
                np.abs((T_fuel_new - T_fuel[aid]) / (np.abs(T_fuel[aid]) + eps))))
            max_dT_inner_rel = max(max_dT_inner_rel, np.max(
                np.abs((T_inner_new - T_inner[aid]) / (np.abs(T_inner[aid]) + eps))))

            T_fuel[aid]    = T_fuel_new
            T_clad[aid]    = T_clad_new
            T_inner[aid]   = T_inner_new
            rho_inner[aid] = rho_inner_new

        # ── bypass pool update ────────────────────────────────
        valid = {aid: r for aid, r in th_results.items()
                 if not is_th_result_degenerate(r, config)}

        if valid:
            total_P = sum(r["P_gen"] for r in valid.values())
            if total_P > 1.0:
                T_outer_th   = sum(r["T_outer_z_K"] * r["P_gen"]
                                   for r in valid.values()) / total_P
                T_pool_mean  = float(np.mean(T_outer_th))
            else:
                T_pool_mean = T_init
        else:
            T_pool_mean = T_init

        T_outer   = np.full(config.n_z, T_pool_mean)
        rho_outer = water_density_iapws(T_outer)

        # ── power convergence ─────────────────────────────────
        sorted_ids     = sorted(valid.keys())
        P_total_new    = sum(th_results[aid]["P_gen"] for aid in valid)
        P_per_assy_new = np.array([th_results[aid]["P_gen"] for aid in sorted_ids])

        dP_total_rel = abs(P_total_new - P_total_prev) / (abs(P_total_prev) + eps)
        n = min(len(P_per_assy_new), len(P_per_assy_prev))
        dP_assy_rel  = np.max(np.abs(P_per_assy_new[:n] - P_per_assy_prev[:n]) /
                              (np.abs(P_per_assy_prev[:n]) + eps))

        P_total_prev    = P_total_new
        P_per_assy_prev = P_per_assy_new

        converged = (
            it >= 5 and
            dP_total_rel < tol_power_total and
            dP_assy_rel  < tol_power_assy
        )

        # ── report ────────────────────────────────────────────
        print(f"P_total         = {P_total_new:.1f} W")
        print(f"dP_total_rel    = {dP_total_rel:.4e}")
        print(f"dP_assy_rel_max = {dP_assy_rel:.4e}")
        print(f"max fuel  dT rel  = {max_dT_fuel_rel:.6e}")
        print(f"max inner dT rel  = {max_dT_inner_rel:.6e}")
        print(f"T_bypass mean     = {T_pool_mean:.4f} K")

        if log_file:
            log_iteration(log_file, it, th_results)

        os.chdir(base_dir)

        if converged:
            print(f"\nConverged in {it + 1} iterations.")
            break

    print("\nFeedback loop finished.")


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run neutronics-TH feedback loop.")
    parser.add_argument("case",              help=f"Case name. Available: {list(cases)}")
    parser.add_argument("--max-iter", "-n",  type=int,   default=10)
    parser.add_argument("--alpha",    "-a",  type=float, default=1.0)
    parser.add_argument("--mesh-dir", "-m",  default="/home/faticen/master_thesis/coupling/meshes")
    args = parser.parse_args()

    if args.case not in cases:
        raise SystemExit(f"Unknown case '{args.case}'. Available: {list(cases)}")

    config    = cases[args.case]
    base_dir  = os.path.dirname(os.path.abspath(__file__))
    log_path  = os.path.join(base_dir, f"case_{args.case}", "iteration_log.txt")
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    with open(log_path, "w", buffering=1) as log_file:
        full_feedback_loop(
            base_dir,
            args.mesh_dir,
            args.case,
            config,
            max_iter=args.max_iter,
            alpha=args.alpha,
            log_file=log_file,
        )