
import openmc
import numpy as np
import glob
import os
from model_builder import build_model
from cases import cases
import re
import pandas as pd


def flatten_columns(df: pd.DataFrame) -> pd.DataFrame:
    new_cols = []
    for col in df.columns:
        if isinstance(col, tuple):
            clean = "_".join([str(c) for c in col if c != ""])
        else:
            clean = str(col)
        new_cols.append(clean)
    df.columns = new_cols
    return df


def enrich_openmc_df_volume_integrated(df: pd.DataFrame, tally: openmc.Tally):
    """
    Add voxel volume + centroid coords.
    Also create volume-INTEGRATED bin values:
      mean_int = mean * voxel_volume
      std_int  = std  * voxel_volume
    (Assuming mean/std are volume-averaged in each mesh bin, which is typical.)
    """
    df = flatten_columns(df)

    mesh_filter = tally.find_filter(openmc.MeshFilter)
    mesh = mesh_filter.mesh
    volumes = mesh.volumes
    centroids = mesh.centroids  # shape (nx, ny, nz, 3) or (nr,nphi,nz,3) for CylindricalMesh

    # Try to find mesh index columns
    # Common names: "mesh 1_x", "mesh 1_y", "mesh 1_z" OR "mesh 1_r", "mesh 1_phi", "mesh 1_z"
    mesh_cols = [c for c in df.columns if "mesh" in c and (c.endswith("_x") or c.endswith("_y") or c.endswith("_z") or c.endswith("_r") or c.endswith("_phi"))]

    # Identify indices
    # Prefer r/phi/z if present, else x/y/z
    r_like = [c for c in df.columns if c.endswith("_r")]
    phi_like = [c for c in df.columns if c.endswith("_phi")]
    x_like = [c for c in df.columns if c.endswith("_x")]
    y_like = [c for c in df.columns if c.endswith("_y")]
    z_like = [c for c in df.columns if c.endswith("_z")]

    if not z_like:
        raise KeyError("Cannot find mesh axial index column (ending with '_z').")

    z_col = z_like[0]

    if r_like and phi_like:
        r_col = r_like[0]
        phi_col = phi_like[0]
        r_idx = df[r_col].values - 1
        phi_idx = df[phi_col].values - 1
        z_idx = df[z_col].values - 1

        df["voxel_volume"] = volumes[r_idx, phi_idx, z_idx]
        df["x_coord"] = centroids[r_idx, phi_idx, z_idx, 0]
        df["y_coord"] = centroids[r_idx, phi_idx, z_idx, 1]
        df["z_coord"] = centroids[r_idx, phi_idx, z_idx, 2]

    elif x_like and y_like:
        x_col = x_like[0]
        y_col = y_like[0]
        x_idx = df[x_col].values - 1
        y_idx = df[y_col].values - 1
        z_idx = df[z_col].values - 1

        df["voxel_volume"] = volumes[x_idx, y_idx, z_idx]
        df["x_coord"] = centroids[x_idx, y_idx, z_idx, 0]
        df["y_coord"] = centroids[x_idx, y_idx, z_idx, 1]
        df["z_coord"] = centroids[x_idx, y_idx, z_idx, 2]

    else:
        raise KeyError("Could not identify mesh indices (need either r/phi/z or x/y/z columns).")

    return df



def build_qprime_pin_for_solver(statepoint_path, name, config):

    sp_file = openmc.StatePoint(statepoint_path)

    # 1) total kappa-fission (eV/source)
    
    Ktot_eV_per_source = float(np.asarray(sp_file.get_tally(name="global_kappa_fission").mean).sum())
    P_core_W = config.P_core_W # <-- set your modeled power (W)

    print("P core [W] =", P_core_W)
    eV_to_J = 1.602e-19
    Ktot_J_per_source = Ktot_eV_per_source * eV_to_J

    f_src_per_s = P_core_W / Ktot_J_per_source  # [source/s]

    print("Ktot [eV/source] =", Ktot_eV_per_source)
    print("f_src_per_s [1/s] =", f_src_per_s)

    def tally_to_axial(sp, name, N_z):
        tal = sp.get_tally(name=name)
        w = np.asarray(tal.mean).reshape(-1)

        if w.size == N_z:
            w_z = w
        else:
            if w.size % N_z != 0:
                raise ValueError(f"{name}: {w.size} bins not divisible by N_z={N_z}")
            Nr = w.size // N_z
            w_z = w.reshape(Nr, N_z).sum(axis=0)

        return np.clip(w_z, 0.0, None)

    pat = re.compile(r"^kappa_fission_tally_cell_(\d+)$")
    items = [(int(pat.match(t.name).group(1)), t.name) for t in sp_file.tallies.values() if pat.match(t.name)]
    items.sort(key=lambda x: -x[0])  # 193, 192, ...

    n_pins = len(items)
    K_pin_z = np.zeros((n_pins, config.n_z))

    for p in range(n_pins):
        K_pin_z[p, :] = tally_to_axial(sp_file, items[p][1], config.n_z)

    # power per axial bin (W) = K[eV/source] * eV_to_J * f[src/s]
    P_bin_W = K_pin_z * eV_to_J * f_src_per_s

    # convert to W/m (linear power in each bin)
    qprime_pin = P_bin_W / config.dz


    # # group pins by assembly
    # assembly_map = {}

    # for pin_id in range(1, config.n_pins_total + 1):

    #     assembly_id = config.pin_to_assembly[pin_id]

    #     assembly_map.setdefault(assembly_id, [])
    #     assembly_map[assembly_id].append(pin_id)

    # qprime_by_assembly = {}

    # for assy_id, pins in assembly_map.items():

    #     pins_sorted = sorted(pins)

    #     arr = np.zeros((len(pins_sorted), config.n_z))

    #     for i, pin in enumerate(pins_sorted):
    #         arr[i, :] = qprime_pin[pin - 1, :]

    #     qprime_by_assembly[assy_id] = arr

    print("qprime_pin shape:", qprime_pin.shape)
    print("Integrated power [W]:", (qprime_pin.sum() * config.dz))

    print("Saving qprime for case:", name)
    print("Working directory:", os.getcwd())

    np.save("qprime_pin.npy", qprime_pin)

    return qprime_pin


def build_qprime_pin_for_solver_by_assembly(statepoint_path, config):

    sp = openmc.StatePoint(statepoint_path)

    # -----------------------------
    # 1) normalization
    # -----------------------------
    Ktot = float(np.asarray(
        sp.get_tally(name="global_kappa_fission").mean
    ).sum())

    eV_to_J = 1.602e-19
    f_src = config.P_core_W / (Ktot * eV_to_J)

    print("f_src_per_s =", f_src)

    # -----------------------------
    # helper
    # -----------------------------
    def tally_to_axial(tally):
        w = np.asarray(tally.mean).reshape(-1)

        if w.size == config.n_z:
            return w

        if w.size % config.n_z != 0:
            raise ValueError("Bad tally shape")

        Nr = w.size // config.n_z
        return w.reshape(Nr, config.n_z).sum(axis=0)

    # -----------------------------
    # 2) extract tallies
    # -----------------------------
    pattern = re.compile(r"kappa_fission_tally_cell_a(\d+)_p(\d+)")

    qprime_by_assembly = {}

    for tally in sp.tallies.values():

        m = pattern.match(tally.name)
        if not m:
            continue

        assy = int(m.group(1))
        pin  = int(m.group(2))

        vals = tally_to_axial(tally)

        # convert to W
        power_bin = vals * eV_to_J * f_src

        # convert to W/cm
        qprime = power_bin / config.dz

        mesh_filter = tally.find_filter(openmc.MeshFilter)
        mesh = mesh_filter.mesh

        x0, y0, _ = mesh.origin

        qprime_by_assembly.setdefault(assy, {
            "values": [],
            "centers": []
        })

        qprime_by_assembly[assy]["values"].append(qprime)
        qprime_by_assembly[assy]["centers"].append((x0, y0))

    # -----------------------------
    # 3) convert to arrays
    # -----------------------------
    # for assy in qprime_by_assembly:

    #     pins = sorted(qprime_by_assembly[assy].keys())

    #     arr = np.array([qprime_by_assembly[assy][p] for p in pins])

    #     qprime_by_assembly[assy] = arr

    for assy in qprime_by_assembly:

        values = np.array(qprime_by_assembly[assy]["values"])
        centers = np.array(qprime_by_assembly[assy]["centers"])

        # sort spatially
        sort_idx = np.lexsort((centers[:,1], centers[:,0]))
        values = values[sort_idx]
        centers = centers[sort_idx]

        qprime_by_assembly[assy] = {
            "values": values,
            "centers": centers
        }

        print(f"Assembly {assy}: {values.shape[0]} pins detected")

    # -----------------------------
    # 4) check total power
    # -----------------------------
    total_power = 0.0
    for data in qprime_by_assembly.values():
        total_power += data["values"].sum() * config.dz

    print("Recovered core power [W]:", total_power)

    return qprime_by_assembly



def build_qprime_pin_for_solver_by_assembly_new(statepoint_path, config):

    import numpy as np
    import re
    import openmc

    sp = openmc.StatePoint(statepoint_path)

    # -----------------------------
    # 1) normalization
    # -----------------------------
    Ktot = float(np.asarray(
        sp.get_tally(name="global_kappa_fission").mean
    ).sum())

    eV_to_J = 1.602e-19
    f_src = config.P_core_W / (Ktot * eV_to_J)

    print("f_src_per_s =", f_src)


    # -----------------------------
    # 2) loop over mesh tallies
    # -----------------------------
    pattern = re.compile(r"kappa_fission_tally_cell_a(\d+)_p(\d+)")

    qprime_by_assembly = {}

    for tally in sp.tallies.values():

        m = pattern.match(tally.name)
        if not m:
            continue

        assy = int(m.group(1))

        # --- dataframe with geometry
        df = tally.get_pandas_dataframe()
        df = enrich_openmc_df_volume_integrated(df, tally)

        # -----------------------------
        # 3) convert to power (IMPORTANT FIX)
        # -----------------------------
        df["power"] = df["mean"] * eV_to_J * f_src   # W per voxel

        # axial linear heat rate (W/m)
        df["qprime"] = df["power"] / config.dz

        # -----------------------------
        # 4) group by pins (robust)
        # -----------------------------
        # avoid float precision issues
        tol = 1e-5
        df["x_round"] = (df["x_coord"] / tol).round().astype(int)
        df["y_round"] = (df["y_coord"] / tol).round().astype(int)

        pin_groups = df.groupby(["x_round", "y_round"])

        pin_profiles = []
        pin_centers  = []

        for (xr, yr), g in pin_groups:

            # real coordinates (mean)
            x = g["x_coord"].mean()
            y = g["y_coord"].mean()

            g = g.sort_values("z_coord")

            qprime = g["qprime"].values

            # sanity check
            if len(qprime) != config.n_z:
                raise ValueError(
                    f"Bad axial size for assembly {assy}: "
                    f"{len(qprime)} != {config.n_z}"
                )

            pin_profiles.append(qprime)
            pin_centers.append((x, y))

        if len(pin_profiles) == 0:
            raise ValueError(f"No pins detected for assembly {assy}")

        arr = np.array(pin_profiles)
        centers = np.array(pin_centers)

        # -----------------------------
        # 5) sort pins (IMPORTANT for stability)
        # -----------------------------
        sort_idx = np.lexsort((centers[:, 1], centers[:, 0]))
        arr = arr[sort_idx]
        centers = centers[sort_idx]

        if assy not in qprime_by_assembly:
            qprime_by_assembly[assy] = {
                "values": [],
                "centers": []
            }

        qprime_by_assembly[assy]["values"].append(arr)
        qprime_by_assembly[assy]["centers"].append(centers)

        print(f"Assembly {assy}: {arr.shape[0]} pins detected")

    for assy in qprime_by_assembly:

        values_list = qprime_by_assembly[assy]["values"]
        centers_list = qprime_by_assembly[assy]["centers"]

        values = np.vstack(values_list)
        centers = np.vstack(centers_list)

        # sort for stability
        sort_idx = np.lexsort((centers[:,1], centers[:,0]))
        values = values[sort_idx]
        centers = centers[sort_idx]

        qprime_by_assembly[assy] = {
            "values": values,
            "centers": centers
        }

        print(f"Assembly {assy}: {values.shape[0]} pins detected")

    # -----------------------------
    # 6) check total power
    # -----------------------------
    total_power = 0.0
    for data in qprime_by_assembly.values():
        total_power += data["values"].sum() * config.dz

    print("Recovered core power [W]:", total_power)

    return qprime_by_assembly