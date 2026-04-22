import os, glob
import numpy as np
import openmc

from model_builder import build_model
from power_extraction import build_qprime_pin_for_solver_by_assembly
from dual_velocity_solver import Dual_Velocity_Assembly_Solver
from cases import cases
from model_reconstruct import reconstruct_model, build_fuel_slices, build_clad_slices, build_coolant_slices, create_axial_planes
import matplotlib.pyplot as plt
from iapws import IAPWS97

def water_density_simple(T_K):
    rho_kg_m3 = 1000.0 * (1.0 - 0.0003 * (np.asarray(T_K) - 293.15))
    return rho_kg_m3 / 1000.0  # g/cm3


def water_density_iapws(T_K, P_MPa=0.1013):
    T = np.asarray(T_K)
    rho = np.zeros_like(T, dtype=float)

    for i, Ti in np.ndenumerate(T):
        state = IAPWS97(T=Ti, P=P_MPa)
        rho[i] = state.rho  # kg/m3

    return rho / 1000.0  # g/cm3

def flatten_cell_dict(cell_dict):
    flat = []
    for assy_id in sorted(cell_dict.keys()):
        flat.extend(cell_dict[assy_id])
    return flat


def get_mesh_path(mesh_dir, config, assembly_type):
    return os.path.join(mesh_dir, config.mesh_map[assembly_type])

def log_iteration(log_file, it, th_results):
    log_file.write(f"\n=== Iteration {it} ===\n")

    for assy_id in sorted(th_results.keys()):
        res = th_results[assy_id]

        T_inner = res["T_inner_z_K"]
        rho_in = res["rho_inner_z"]

        P_gen = res["P_gen"]
        P_in = res["P_to_inner"]
        P_out = res["P_to_outer"]
        P_mcp = res["P_mcpDT_in"]

        mdot_in = res["mdot_in"]
        mdot_out = res["mdot_out"]


        # frac_in = P_in / P_gen if P_gen > 0 else 0.0
        # frac_out = P_out / P_gen if P_gen > 0 else 0.0

        log_file.write(
            f"Assembly {assy_id}: "
            f"T_in={T_inner[0]:.2f}, "
            f"T_out={T_inner[-1]:.2f}, "
            f"rho_in={rho_in[0]:.5f}, "
            f"rho_out={rho_in[-1]:.5f}, "
            f"P_gen={P_gen:.1f}, "
            f"P_in={P_in:.1f}, "
            f"P_out={P_out:.1f}, "
            f"P_mcp={P_mcp:.1f}, "
            # f"f_in={frac_in:.2f}, "
            # f"f_out={frac_out:.2f}\n"
            f"mdot_in={mdot_in:.4f}, "
            f"mdot_out={mdot_out:.4f}"
        )

    log_file.flush()

def run_openmc_in_dir(model, workdir, PlotColors):

    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    # export model first
    model.export_to_xml()
    # --- DEBUG PLOT
    plot = openmc.Plot()
    plot.color_by = 'material'
    plot.basis = 'xz'
    plot.origin = (0, -17.15, 0)
    plot.width = (32, 75)
    plot.pixels = (1500,1500)
    plot.colors = PlotColors

    plot_1 = openmc.Plot()
    plot_1.color_by = 'cell'
    plot_1.basis = 'xy'
    plot_1.origin = (0, -18, 0)
    plot_1.width = (32, 32)
    plot_1.pixels = (1500,1500)

    plot_2 = openmc.Plot()
    plot_2.color_by = 'material'
    plot_2.basis = 'xz'
    plot_2.origin = (0, -15.45, 0)
    plot_2.width = (32, 75)
    plot_2.pixels = (1500,1500)
    plot_2.colors = PlotColors

    plot_3 = openmc.Plot()
    plot_3.color_by = 'material'
    plot_3.basis = 'xz'
    plot_3.origin = (0, -18.85, 0)
    plot_3.width = (32, 75)
    plot_3.pixels = (1500,1500)
    plot_3.colors = PlotColors


    plot_4 = openmc.Plot()
    plot_4.color_by = 'material'
    plot_4.basis = 'xz'
    plot_4.origin = (0, 0, 0)
    plot_4.width = (32, 75)
    plot_4.pixels = (1500,1500)
    plot_4.colors = PlotColors


    plots = openmc.Plots([plot, plot_1, plot_2, plot_3, plot_4])
    plots.export_to_xml()

    openmc.plot_geometry()

    print(model.geometry.bounding_box)

    
    openmc.run()

    sp_path = sorted(glob.glob("statepoint.*.h5"))[-1]
    return os.path.join(workdir, sp_path)



def run_th_case(config, mesh_path, qprime_data, case_name, outdir):
    """
    Updated to use structured qprime_data with automatic center matching.
    
    Args:
        qprime_data: dict with 'values' (n_pins, N_z) and 'centers' (n_pins, 2)
    """
    solver = Dual_Velocity_Assembly_Solver(config, mesh_path)
    qprime_pin = qprime_data['values']
    rod_centers = qprime_data['centers'] 
    solver.assign_qprime_by_centers(qprime_pin, rod_centers)  
    solver.set_power_profile(qprime_pin)
    solver.run_old(case_name=case_name, out_dir=outdir)

    T_inner_z_K = np.asarray(solver.plot_T_bulk_in) + 273.15
    T_outer_z_K = np.asarray(solver.plot_T_bulk_out) + 273.15
    T_fuel_z_K  = np.asarray(solver.plot_T_fuel) + 273.15
    T_clad_z_K  = 0.5 * (T_fuel_z_K + T_inner_z_K)

    rho_inner_z = water_density_iapws(T_inner_z_K)
    rho_outer_z = water_density_iapws(T_outer_z_K)

    mdot_in = solver.plot_mdot_in
    mdot_out = solver.plot_mdot_out

    P_gen = np.sum(solver.qprime_pin) * solver.dz
    P_to_inner = np.sum(solver.plot_heat_in) * solver.dz
    P_to_outer = np.sum(solver.plot_heat_out) * solver.dz

    dT = np.diff(solver.plot_T_bulk_in, prepend=solver.plot_T_bulk_in[0])
    P_mcpDT_in = np.sum(solver.plot_mdot_in * solver.plot_cp_in * dT)


    print(f"P_gen           = {P_gen:.3f} W")
    print(f"P_to_inner      = {P_to_inner:.3f} W")
    print(f"P_to_outer      = {P_to_outer:.3f} W")
    print(f"P_inner+outer   = {P_to_inner + P_to_outer:.3f} W")
    print(f"m_dot*cp*dT in  = {P_mcpDT_in:.3f} W")
    

    return {
    "qprime_pin": qprime_pin,
    "T_inner_z_K": T_inner_z_K,
    "T_outer_z_K": T_outer_z_K,
    "T_fuel_z_K": T_fuel_z_K,
    "T_clad_z_K": T_clad_z_K,
    "rho_inner_z": rho_inner_z,
    "rho_outer_z": rho_outer_z,
    "P_gen": P_gen,
    "P_to_inner": P_to_inner,
    "P_to_outer": P_to_outer,
    "P_mcpDT_in": P_mcpDT_in,
    "mdot_in": mdot_in,
    "mdot_out": mdot_out
    }

def full_feedback_loop(base_dir, mesh_path, case_name, config, n_iter=6):
    os.chdir(base_dir)

    model, geometry_data, outer_coolant_cells, PlotColors = build_model(config)

    universe_0 = model.geometry.root_universe
    print('Geomtery unchanged:', universe_0)
    
    z_planes = create_axial_planes(config.z_min, config.z_max, config.n_z)

    fuel_cells_by_assembly = geometry_data["fuel_cells"]
    clad_cells_by_assembly = geometry_data["clad_cells"]
    inner_coolant_cells_by_assembly = geometry_data["inner_coolant_cells"]
    rod_mapping = geometry_data["rod_mapping"]
    assembly_types = geometry_data["assembly_types_by_id_openmc"]
    detector_rod_centers = geometry_data["detector_rod_centers"]

    water_below_cells_by_assembly = geometry_data["water_below_assembly_cells"]
    water_above_cells_by_assembly = geometry_data["water_above_assembly_cells"]
    top_grid_cells_by_assembly = geometry_data["top_grid_cells"]
    bottom_grid_cells_by_assembly = geometry_data["bottom_grid_cells"]
    can_cells_by_assembly = geometry_data["can_cells"]
    can_bottom_walls_cells_by_assembly = geometry_data["can_bottom_walls_cells"]

    fuel_cells = flatten_cell_dict(fuel_cells_by_assembly)
    clad_cells = flatten_cell_dict(clad_cells_by_assembly)
    inner_coolant_cells = flatten_cell_dict(inner_coolant_cells_by_assembly)

    fuel_layers_by_assembly = {}
    clad_layers_by_assembly = {}
    inner_coolant_layers_by_assembly = {}
    inner_coolant_materials = {}

    for assy_id, cells in fuel_cells_by_assembly.items():

        fuel_layers_by_assembly[assy_id] = build_fuel_slices(
            cells,
            universe_0,
            z_planes,
            config.n_z
        )

    for assy_id, cells in clad_cells_by_assembly.items():

        layers, lower, upper = build_clad_slices(
            cells,
            universe_0,
            z_planes,
            config.n_z
        )

        clad_layers_by_assembly[assy_id] = {
            "layers": layers,
            "lower": lower,
            "upper": upper
        }

    for assy_id, cells in inner_coolant_cells_by_assembly.items():    

        layers, lower_cell, upper_cell = build_coolant_slices(
            cells,
            universe_0,
            z_planes,
            config.n_z
        )

        inner_coolant_layers_by_assembly[assy_id] = {
            "layers": layers,
            "lower": lower_cell,
            "upper": upper_cell
        }

        inner_coolant_materials[assy_id] = []

        # take representative material
        base_mat = cells[0].fill

        for k in range(config.n_z):
            m = base_mat.clone()
            m.name = f"coolant_a{assy_id}_z{k}"
            inner_coolant_materials[assy_id].append(m)



    for assy_id, data in inner_coolant_layers_by_assembly.items():

        layers = data["layers"]

        for k in range(config.n_z):
            for c in layers[k]:
                c.fill = inner_coolant_materials[assy_id][k]



    inner_coolant_boundary_materials = {}

    for assy_id, cells in inner_coolant_cells_by_assembly.items():

        base_mat = cells[0].fill

        lower_mat = base_mat.clone()
        lower_mat.name = f"coolant_a{assy_id}_lower"

        upper_mat = base_mat.clone()
        upper_mat.name = f"coolant_a{assy_id}_upper"

        inner_coolant_boundary_materials[assy_id] = {
            "lower": lower_mat,
            "upper": upper_mat
        }


    for assy_id, data in inner_coolant_layers_by_assembly.items():

        data["lower"].fill = inner_coolant_boundary_materials[assy_id]["lower"]
        data["upper"].fill = inner_coolant_boundary_materials[assy_id]["upper"]



    all_inner_coolant_materials = []

    for mats in inner_coolant_materials.values():
        all_inner_coolant_materials.extend(mats)

    print(all_inner_coolant_materials)

    model.materials.extend(all_inner_coolant_materials)

    all_boundary_materials = []

    for mats in inner_coolant_boundary_materials.values():
        all_boundary_materials.extend(mats.values())


    print(all_boundary_materials)

    model.materials.extend(all_boundary_materials)
 

   
    outer_coolant_layers, outer_lower, outer_upper = build_coolant_slices(
        outer_coolant_cells,
        universe_0,
        z_planes,
        config.n_z
    )

    outer_coolant_data = {
        "layers": outer_coolant_layers,
        "lower": outer_lower,
        "upper": outer_upper
    }


    # After inner coolant materials setup, add:
    outer_coolant_materials = {
        "layers": {},
        "lower": None,
        "upper": None
    }

    base_outer_mat = outer_coolant_cells[0].fill

    # --- axial layers
    for k in range(config.n_z):
        m = base_outer_mat.clone()
        m.name = f"outer_coolant_z{k}"
        outer_coolant_materials["layers"][k] = m

    # --- lower region
    m_lower = base_outer_mat.clone()
    m_lower.name = "outer_coolant_lower"
    outer_coolant_materials["lower"] = m_lower

    # --- upper region
    m_upper = base_outer_mat.clone()
    m_upper.name = "outer_coolant_upper"
    outer_coolant_materials["upper"] = m_upper

    # layers
    for k in range(config.n_z):
        for c in outer_coolant_layers[k]:
            c.fill = outer_coolant_materials["layers"][k]

    # lower
    outer_lower.fill = outer_coolant_materials["lower"]

    # upper
    outer_upper.fill = outer_coolant_materials["upper"]
        
    # layers (list)
    model.materials.extend(outer_coolant_materials["layers"])

    # single materials
    model.materials.append(outer_coolant_materials["lower"])
    model.materials.append(outer_coolant_materials["upper"]) 



    # --- water below
    water_above_materials = {}
    water_below_materials = {}

    for assy_id, cells in water_above_cells_by_assembly.items():
        water_above_materials[assy_id] = []
        for cell in cells:
            mat = cell.fill.clone()
            cell.fill = mat
            water_above_materials[assy_id].append(mat)

    for assy_id, cells in water_below_cells_by_assembly.items():
        water_below_materials[assy_id] = []
        for cell in cells:
            mat = cell.fill.clone()
            cell.fill = mat
            water_below_materials[assy_id].append(mat)

    all_above = [m for mats in water_above_materials.values() for m in mats]
    all_below = [m for mats in water_below_materials.values() for m in mats]

    model.materials.extend(all_above)
    model.materials.extend(all_below)

    # print('Geomtery divided:', universe_0)


    # Initial uniform guesses
    T_init = getattr(config, "T_inlet_K", 295.65)
    rho_innit =  water_density_iapws(T_init)

    T_fuel = {}
    T_clad = {}
    T_inner = {}
    rho_inner = {}

    for assy_id in fuel_cells_by_assembly:

        T_fuel[assy_id] = np.full(config.n_z, T_init)
        T_clad[assy_id] = np.full(config.n_z, T_init)
        T_inner[assy_id] = np.full(config.n_z, T_init)
        rho_inner[assy_id] = np.full(config.n_z, water_density_iapws(T_init))

    # Add to initial state dicts
    T_outer = np.full(config.n_z, T_init)
    rho_outer = np.full(config.n_z, water_density_simple(T_init))

    max_iter = 30
    tol_T = 1e-2
    tol_rho = 1e-3
    alpha = 0.2


    # Before the iteration loop, make sure:
    th = {
        "T_fuel": {assy_id: np.zeros(config.n_z) for assy_id in fuel_cells_by_assembly},
        "T_clad": {assy_id: np.zeros(config.n_z) for assy_id in fuel_cells_by_assembly},
        "T_inner": {assy_id: np.zeros(config.n_z) for assy_id in fuel_cells_by_assembly},
        "rho_inner": {assy_id: np.zeros(config.n_z) for assy_id in fuel_cells_by_assembly}
    }

    for it in range(max_iter):
        print(f"\n================ ITERATION {it} ================")

        for assy_id in fuel_layers_by_assembly:

            for k in range(config.n_z):

                for c in fuel_layers_by_assembly[assy_id][k]:
                    c.temperature = float(T_fuel[assy_id][k])


        for assy_id, data in clad_layers_by_assembly.items():

            # active region
            for k, cells in enumerate(data["layers"]):
                for c in cells:
                    c.temperature = T_clad[assy_id][k]

            # lower region (no TH → assume inlet or constant)
            data["lower"].temperature = T_clad[assy_id][0]

            # upper region
            data["upper"].temperature = T_clad[assy_id][-1]


        # INNER COOLANT
        for assy_id, data in inner_coolant_layers_by_assembly.items():

            for k, cells in enumerate(data["layers"]):
                for c in cells:
                    c.temperature = float(T_inner[assy_id][k])

                inner_coolant_materials[assy_id][k].set_density(
                    "g/cm3",
                    float(rho_inner[assy_id][k])
                )

            # lower
            data["lower"].temperature = float(T_init)

            inner_coolant_boundary_materials[assy_id]["lower"].set_density(
                "g/cm3",
                rho_innit)
            

            # upper
            data["upper"].temperature = float(T_inner[assy_id][-1])

            inner_coolant_boundary_materials[assy_id]["upper"].set_density(
                "g/cm3",
                float(rho_inner[assy_id][-1])
            )


            T_top = float(T_inner[assy_id][-1])
            rho_top = float(rho_inner[assy_id][-1])

            T_bottom = float(T_inner[assy_id][0])
            rho_bottom = float(rho_inner[assy_id][0])

            # T_can_z = 0.5 * (T_inner_z + T_outer_z)
            # T_can_avg = float(np.mean(T_can_z))

            T_can_bottom = float(T_inner[assy_id][0])

            # --- water above
            for cell in water_above_cells_by_assembly[assy_id]:
                cell.temperature = T_top
                cell.fill.set_density("g/cm3", rho_top)

            # --- water below
            for cell in water_below_cells_by_assembly[assy_id]:
                cell.temperature = T_bottom   # or T_bottom
                cell.fill.set_density("g/cm3", rho_bottom)

            # --- top grid (homogenized al+water)
            for cell in top_grid_cells_by_assembly[assy_id]:
                cell.temperature = T_top
                # only if density feedback makes sense for this homogenized material:
                # update_homogenized_grid_density(cell, rho_top)

            # --- bottom grid
            for cell in bottom_grid_cells_by_assembly[assy_id]:
                cell.temperature = T_bottom
                # update_homogenized_grid_density(cell, rho_bottom)

            

        # OUTER COOLANT
        for k in range(config.n_z):
            for c in outer_coolant_layers[k]:
                c.temperature = float(T_outer[k])

            outer_coolant_materials["layers"][k].set_density(
                "g/cm3",
                float(rho_outer[k])
            )

        # lower
        outer_coolant_data["lower"].temperature = float(T_outer[0])

        outer_coolant_materials["lower"].set_density(
            "g/cm3",
            float(rho_outer[0])
        )

        # upper
        outer_coolant_data["upper"].temperature = float(T_outer[-1])

        outer_coolant_materials["upper"].set_density(
            "g/cm3",
            float(rho_outer[-1])
        )
      
        model.materials = openmc.Materials(
            model.geometry.get_all_materials().values()
        )
        # print(model.materials)


        if it == 0:
            particles = 10000
            batches = 200
            inactive = 50
        elif it <= 3:
            particles = 30000
            batches = 300
            inactive = 80
        elif it <= 6:
            particles = 60000
            batches = 600
            inactive = 100
        elif it <= 10:
            particles = 80000
            batches = 800
            inactive = 100
        else:
            particles = 100000
            batches = 1000
            inactive = 200

        model.settings.particles = particles
        model.settings.batches = batches
        model.settings.inactive = inactive

        print(f"OpenMC stats: particles={particles}, batches={batches}")


        # Run OpenMC
        omc_dir = os.path.join(base_dir, f"case_{case_name}", f"it_{it:02d}", "openmc")
        sp_path = run_openmc_in_dir(model, omc_dir, PlotColors=PlotColors)

        qprime_by_assembly = build_qprime_pin_for_solver_by_assembly(sp_path, config)

        qprime_data = {}
        for assy_id, qprime_array in qprime_by_assembly.items():
            # Get rod centers for this assembly from rod_mapping
            assembly_rods = rod_mapping[assy_id]
            centers = np.array([rod['center'] for rod in assembly_rods])
            
            # Store structured data for TH solver
            qprime_data[assy_id] = {
                'values': qprime_array,  # q' matrix from tallies (W/m per pin)
                'centers': centers       # (N_pins, 2) array of (x,y) rod centers
            }

        # Run TH
        th_dir = os.path.join(base_dir, f"case_{case_name}", f"it_{it:02d}", "th")
        os.makedirs(th_dir, exist_ok=True)

        th_results = {}

        for assy_id, assy_type in assembly_types.items():

            print(f"Running TH for assembly {assy_id} ({assy_type})")

            mesh_file = get_mesh_path(mesh_path, config, assy_type)

            qprime_data_assy = qprime_data[assy_id]     # Power profile (W/m)

            th_results[assy_id] = run_th_case(
                config,
                mesh_file,
                qprime_data_assy,   #  IMPORTANT CHANGE
                case_name=f"{case_name}_it{it:02d}_assy{assy_id}",
                outdir=os.path.join(th_dir, f"assy_{assy_id}")
            )


        tol_fuel_abs = 0.05      # K
        tol_fuel_rel = 1e-4

        tol_innerT_abs = 0.05    # K
        tol_innerT_rel = 1e-4

        tol_rho_abs = 1e-4       # g/cm3
        tol_rho_rel = 1e-4

        tol_outerT_abs = 0.05    # K
        tol_outerT_rel = 1e-4

        tol_outerrho_abs = 1e-4  # g/cm3
        tol_outerrho_rel = 1e-4

        eps = 1e-8

        max_dT_fuel_abs = 0.0
        max_dT_fuel_rel = 0.0

        max_dT_inner_abs = 0.0
        max_dT_inner_rel = 0.0

        max_drho_inner_abs = 0.0
        max_drho_inner_rel = 0.0

        # -----------------------------
        # update assembly-wise fields
        # -----------------------------
        for assy_id, res in th_results.items():
            T_fuel_th = res["T_fuel_z_K"]
            T_clad_th = res["T_clad_z_K"]
            T_inner_th = res["T_inner_z_K"]
            rho_inner_th = res["rho_inner_z"]

            T_fuel_new = alpha * T_fuel_th + (1.0 - alpha) * T_fuel[assy_id]
            T_clad_new = alpha * T_clad_th + (1.0 - alpha) * T_clad[assy_id]
            T_inner_new = alpha * T_inner_th + (1.0 - alpha) * T_inner[assy_id]
            rho_inner_new = alpha * rho_inner_th + (1.0 - alpha) * rho_inner[assy_id]

            dT_fuel_abs = np.max(np.abs(T_fuel_new - T_fuel[assy_id]))
            dT_fuel_rel = np.max(np.abs((T_fuel_new - T_fuel[assy_id]) / (np.abs(T_fuel[assy_id]) + eps)))

            dT_inner_abs = np.max(np.abs(T_inner_new - T_inner[assy_id]))
            dT_inner_rel = np.max(np.abs((T_inner_new - T_inner[assy_id]) / (np.abs(T_inner[assy_id]) + eps)))

            drho_inner_abs = np.max(np.abs(rho_inner_new - rho_inner[assy_id]))
            drho_inner_rel = np.max(np.abs((rho_inner_new - rho_inner[assy_id]) / (np.abs(rho_inner[assy_id]) + eps)))

            max_dT_fuel_abs = max(max_dT_fuel_abs, dT_fuel_abs)
            max_dT_fuel_rel = max(max_dT_fuel_rel, dT_fuel_rel)

            max_dT_inner_abs = max(max_dT_inner_abs, dT_inner_abs)
            max_dT_inner_rel = max(max_dT_inner_rel, dT_inner_rel)

            max_drho_inner_abs = max(max_drho_inner_abs, drho_inner_abs)
            max_drho_inner_rel = max(max_drho_inner_rel, drho_inner_rel)

            T_fuel[assy_id] = T_fuel_new
            T_clad[assy_id] = T_clad_new
            T_inner[assy_id] = T_inner_new
            rho_inner[assy_id] = rho_inner_new

        # -----------------------------
        # outer coolant: shared field
        # do this ONCE, after assembly loop
        # -----------------------------
        # outer_T_all = np.array([res["T_outer_z_K"] for res in th_results.values()])
        # T_outer_th = np.mean(outer_T_all, axis=0)

        # rho_outer_th = water_density_iapws(T_outer_th)

        # T_outer_new = alpha * T_outer_th + (1.0 - alpha) * T_outer
        # rho_outer_new = alpha * rho_outer_th + (1.0 - alpha) * rho_outer

        # max_dT_outer_abs = np.max(np.abs(T_outer_new - T_outer))
        # max_dT_outer_rel = np.max(np.abs((T_outer_new - T_outer) / (np.abs(T_outer) + eps)))

        # max_drho_outer_abs = np.max(np.abs(rho_outer_new - rho_outer))
        # max_drho_outer_rel = np.max(np.abs((rho_outer_new - rho_outer) / (np.abs(rho_outer) + eps)))

        # T_outer = T_outer_new
        # rho_outer = rho_outer_new

        T_outer_ref = np.full_like(T_outer, T_init)   # or config.T_inlet_K
        rho_outer_ref = water_density_iapws(T_outer_ref)

        T_outer_new = T_outer_ref
        rho_outer_new = rho_outer_ref

        max_dT_outer_abs = np.max(np.abs(T_outer_new - T_outer))
        max_dT_outer_rel = np.max(np.abs((T_outer_new - T_outer) / (np.abs(T_outer) + eps)))

        max_drho_outer_abs = np.max(np.abs(rho_outer_new - rho_outer))
        max_drho_outer_rel = np.max(np.abs((rho_outer_new - rho_outer) / (np.abs(rho_outer) + eps)))

        T_outer = T_outer_new
        rho_outer = rho_outer_new

        print(f"T_outer:", T_outer)


        print(f"max fuel dT abs   = {max_dT_fuel_abs:.6e} K")
        print(f"max fuel dT rel   = {max_dT_fuel_rel:.6e}")

        print(f"max inner dT abs  = {max_dT_inner_abs:.6e} K")
        print(f"max inner dT rel  = {max_dT_inner_rel:.6e}")

        print(f"max inner drho abs= {max_drho_inner_abs:.6e} g/cm3")
        print(f"max inner drho rel= {max_drho_inner_rel:.6e}")

        print(f"max outer dT abs  = {max_dT_outer_abs:.6e} K")
        print(f"max outer dT rel  = {max_dT_outer_rel:.6e}")

        print(f"max outer drho abs= {max_drho_outer_abs:.6e} g/cm3")
        print(f"max outer drho rel= {max_drho_outer_rel:.6e}")

        converged = (
            (max_dT_fuel_abs   < tol_fuel_abs   or max_dT_fuel_rel   < tol_fuel_rel) and
            (max_dT_inner_abs  < tol_innerT_abs or max_dT_inner_rel  < tol_innerT_rel) and
            (max_drho_inner_abs < tol_rho_abs   or max_drho_inner_rel < tol_rho_rel) #and
            # (max_dT_outer_abs  < tol_outerT_abs or max_dT_outer_rel  < tol_outerT_rel) and
            # (max_drho_outer_abs < tol_outerrho_abs or max_drho_outer_rel < tol_outerrho_rel)
        )
        
        log_iteration(log_file, it, th_results)
        os.chdir(base_dir)

        if converged:
            print(f"\nConverged in {it+1} iterations.")
            break

        # for assy_id, res in th_results.items():
        #     th["T_fuel"][assy_id] = res["T_fuel_z_K"]
        #     th["T_clad"][assy_id] = res["T_clad_z_K"]
        #     th["T_inner"][assy_id] = res["T_inner_z_K"]
        #     th["rho_inner"][assy_id] = res["rho_inner_z"]

        #     T_fuel_new = alpha * th["T_fuel"][assy_id] + (1.0 - alpha) * T_fuel[assy_id]
        #     T_clad_new = alpha * th["T_clad"][assy_id] + (1.0 - alpha) * T_clad[assy_id]
        #     T_inner_new = alpha * th["T_inner"][assy_id] + (1.0 - alpha) * T_inner[assy_id]
        #     rho_inner_new = alpha * th["rho_inner"][assy_id] + (1.0 - alpha) * rho_inner[assy_id]

        #     max_dT = max(
        #         max_dT,
        #         np.max(np.abs(T_fuel_new - T_fuel[assy_id])),
        #         np.max(np.abs(T_clad_new - T_clad[assy_id])),
        #         np.max(np.abs(T_inner_new - T_inner[assy_id]))
        #     )

        #     max_drho = max(
        #         max_drho,
        #         np.max(np.abs(rho_inner_new - rho_inner[assy_id]))
        #     )

        #     eps = 1e-8
        #     dT_abs = np.max(np.abs(T_inner_new - T_inner[assy_id]))            
        #     dT_rel = np.max(np.abs((T_inner_new - T_inner[assy_id]) / (T_inner[assy_id] + eps)))

        #     drho_abs = np.max(np.abs(rho_inner_new - rho_inner[assy_id]))
        #     drho_rel = np.max(np.abs((rho_inner_new - rho_inner[assy_id]) / (rho_inner[assy_id] + eps)))

        #     # converged = (
        #     #     (dT_abs < 0.1 or dT_rel < 1e-4) and
        #     #     (drho_abs < 1e-3 or drho_rel < 1e-3)
        #     # )


        #     T_fuel[assy_id] = T_fuel_new
        #     T_clad[assy_id] = T_clad_new
        #     T_inner[assy_id] = T_inner_new
        #     rho_inner[assy_id] = rho_inner_new


        #     # ---- outer coolant: one shared axial profile for whole bypass ----
        #     outer_T_all = np.array([res["T_outer_z_K"] for res in th_results.values()])   # shape (n_assy, n_z)
        #     outer_rho_all = np.array([res["rho_outer_z"] for res in th_results.values()]) # shape (n_assy, n_z)

        #     T_outer_th = np.mean(outer_T_all, axis=0)
        #     rho_outer_th = water_density_simple(T_outer_th)

        #     T_outer_new = alpha * T_outer_th + (1.0 - alpha) * T_outer
        #     rho_outer_new = alpha * rho_outer_th + (1.0 - alpha) * rho_outer

        #     # max_dT = max(max_dT, np.max(np.abs(T_outer_new - T_outer)))
        #     # max_drho = max(max_drho, np.max(np.abs(rho_outer_new - rho_outer)))

        #     T_outer = T_outer_new
        #     rho_outer = rho_outer_new

        #     # T_outer_new = alpha * np.mean([res["T_outer_z_K"] for res in th_results.values()]) + \
        #     #       (1.0 - alpha) * T_outer_assembly_avg  # Track previous avg
    
        #     # rho_outer_new = alpha * np.mean([res["rho_outer_z"] for res in th_results.values()]) + \
        #     #                 (1.0 - alpha) * rho_outer_assembly_avg
            
        #     # # Update axial slices
        #     # for k in range(config.n_z):
        #     #     T_outer[k] = T_outer_new  # Uniform across assemblies per slice
        #     #     rho_outer[k] = rho_outer_new
                
        #     #     max_dT = max(max_dT, abs(T_outer_new - T_outer_assembly_avg))
        #     #     max_drho = max(max_drho, abs(rho_outer_new - rho_outer_assembly_avg))
            
        #     # # Store for next iteration
        #     # T_outer_assembly_avg = T_outer_new
        #     # rho_outer_assembly_avg = rho_outer_new

        #     # # Relax feedback fields
        #     # T_fuel_z  = alpha * th["T_fuel_z_K"]  + (1.0 - alpha) * T_fuel_z
        #     # T_clad_z  = alpha * th["T_clad_z_K"]  + (1.0 - alpha) * T_clad_z
        #     # T_inner_z = alpha * th["T_inner_z_K"] + (1.0 - alpha) * T_inner_z
        #     # T_outer_z = alpha * th["T_outer_z_K"] + (1.0 - alpha) * T_outer_z

        #     # rho_inner_z = alpha * th["rho_inner_z"] + (1.0 - alpha) * rho_inner_z
        #     # rho_outer_z = alpha * th["rho_outer_z"] + (1.0 - alpha) * rho_outer_z

        #     # print("Mean fuel T [K]:", float(np.mean(T_fuel_z)))
        #     # print("Mean inner coolant T [K]:", float(np.mean(T_inner_z)))
        #     # print("Mean outer coolant T [K]:", float(np.mean(T_outer_z)))
        #     # print("Mean inner rho [g/cc]:", float(np.mean(rho_inner_z)))
        #     # print("Mean outer rho [g/cc]:", float(np.mean(rho_outer_z)))
        
        # print(f"max dT = {max_dT:.6e} K")
        # print(f"max drho = {max_drho:.6e} g/cm3")


        # log_iteration(log_file, it, th_results)
        # os.chdir(base_dir)


        # if (dT_abs < 0.1 or dT_rel < 1e-4) and (drho_abs < 1e-3 or drho_rel < 1e-3):
        #     print(f"\nConverged in {it+1} iterations.")
        #     break

    print("\nFeedback loop finished.")



if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    mesh_path = "/home/f8/praksaf8/enidafatic/BME_OpenMC_original/coupling/meshes"

    case_name = "mid_power"
    config = cases[case_name]


    log_path = os.path.join("/home/f8/praksaf8/enidafatic/BME_OpenMC_original/coupling/", f"case_{case_name}", "iteration_log.txt")

    log_file = open(log_path, "w", buffering=1)  # line-buffered


    full_feedback_loop(
        base_dir,
        mesh_path,
        case_name,
        config,
        n_iter=3
    )

