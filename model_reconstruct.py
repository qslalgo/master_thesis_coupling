import openmc 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cases
from model_builder import build_model
import os

def split_cells_into_axial_layers(cells, z_edges, name_prefix):
    """
    Returns layers[k] = list of new OpenMC cells for axial layer k.
    Each original cell is replaced by N_z sliced cells (same fill), with extra z-region.
    """
    N_z = len(z_edges) - 1

    # Create z-planes
    z_planes = []
    for k, z in enumerate(z_edges):
        # boundary_type only on the ends if you want; otherwise leave None
        z_planes.append(openmc.ZPlane(z0=float(z), name=f"{name_prefix}_z{k}"))

    layers = [[] for _ in range(N_z)]
    sliced_cells_all = []

    for base_cell in cells:
        base_region = base_cell.region
        base_fill   = base_cell.fill
        base_temp   = getattr(base_cell, "temperature", None)

        # IMPORTANT: remove base cell from geometry later; we create replacements
        for k in range(N_z):
            zreg = (+z_planes[k] & -z_planes[k+1])

            new_cell = openmc.Cell(
                name=f"{base_cell.name}_{name_prefix}_k{k:03d}",
                fill=base_fill,
                region=base_region & zreg
            )
            if base_temp is not None:
                new_cell.temperature = base_temp

            layers[k].append(new_cell)
            sliced_cells_all.append(new_cell)

    return layers, sliced_cells_all




def reconstruct_model(model, universe_0,
                      fuel_cells, clad_cells,
                      inner_coolant_cells,
                      outer_coolant_cells,
                      config):

    z_min = -25 
    z_max =  25
    N_z   = config.n_z

    z_edges = np.linspace(z_min, z_max, N_z + 1)

    # --- Split cells
    fuel_layers, fuel_sliced = split_cells_into_axial_layers(
        fuel_cells, z_edges, "fuel"
    )

    # clad_layers, clad_sliced = split_cells_into_axial_layers(
    #     clad_cells, z_edges, "clad"
    # )

    #inner_layers, inner_sliced = split_cells_into_axial_layers(
    #    inner_coolant_cells, z_edges, "coolant_inner"
    #)

    #outer_layers, outer_sliced = split_cells_into_axial_layers(
    #    outer_coolant_cells, z_edges, "coolant_outer"
    #)

    # --- Remove base cells from universe
    for cell in fuel_cells:
        if cell.id in universe_0.cells:
            del universe_0.cells[cell.id]

    # --- Add sliced cells
    for cell in fuel_sliced:
        universe_0.add_cell(cell)

    return (
        fuel_layers,
        # clad_layers,
        # inner_layers,
        #outer_layers
    )




def create_axial_planes(zmin, zmax, n_z):
    z_edges = np.linspace(zmin, zmax, n_z + 1)

    planes = []
    for z in z_edges:
        planes.append(
            openmc.ZPlane(
                z0=z
            )
        )

    return planes

def build_fuel_slices(
        fuel_cells,
        root_uni,
        z_planes,
        n_z,
):

    fuel_layers = [[] for _ in range(n_z)]

    for base_fuel_cell in fuel_cells:

        base_region = base_fuel_cell.region
        fuel_material = base_fuel_cell.fill
        base_name = base_fuel_cell.name

        slice_regions = []

        # --- create slices
        for k in range(n_z):

            region = base_region & +z_planes[k] & -z_planes[k+1]
            slice_regions.append(region)

            cell = openmc.Cell(
                fill=fuel_material,
                region=region,
                name=f"{base_name}_z{k}"
            )

            root_uni.add_cell(cell)
            fuel_layers[k].append(cell)

        # --- union of slices
        sliced_union = slice_regions[0]
        for r in slice_regions[1:]:
            sliced_union |= r

        # --- remainder
        remainder_region = base_region & ~sliced_union

        remainder_cell = openmc.Cell(
            fill=fuel_material,
            region=remainder_region,
            name=f"{base_name}_remainder"
        )

        # root_uni.add_cell(remainder_cell)

        if base_fuel_cell.id in root_uni.cells:
            del root_uni.cells[base_fuel_cell.id]

    return fuel_layers


def build_clad_slices(
    clad_cells,
    universe,
    z_planes,
    n_z
):

    clad_layers = [[] for _ in range(n_z)]

    for base_clad_cell in clad_cells:

        base_region = base_clad_cell.region
        base_mat = base_clad_cell.fill
        base_name = base_clad_cell.name

        slice_regions = []

        for k in range(n_z):

            region = base_region & +z_planes[k] & -z_planes[k+1]
            slice_regions.append(region)

            c = openmc.Cell(
                fill=base_mat,
                region=region,
                name=f"{base_name}_z{k}"
            )

            universe.add_cell(c)
            clad_layers[k].append(c)

            # remainder for caps / leftover geometry
            sliced_union = slice_regions[0]
            for r in slice_regions[1:]:
                sliced_union |= r

            remainder_region = base_region & ~sliced_union

            remainder_cell = openmc.Cell(
                fill=base_mat,
                region=remainder_region,
                name=f"{base_name}_remainder"
            )

            # universe.add_cell(remainder_cell)

            lower_region = base_region & -z_planes[0]

            lower_cell = openmc.Cell(
                fill=base_mat,
                region=lower_region,
                name=f"{base_name}_lower"
            )
            universe.add_cell(lower_cell)

            # --- upper coolant (above heated region)
            upper_region = base_region & +z_planes[n_z]

            upper_cell = openmc.Cell(
                fill=base_mat,
                region=upper_region,
                name=f"{base_name}_upper"
            )
            universe.add_cell(upper_cell)

            # remove original clad cell
            if base_clad_cell.id in universe.cells:
                del universe.cells[base_clad_cell.id]

    return clad_layers, lower_cell, upper_cell



def build_coolant_slices(
    coolant_cells,
    universe,
    z_planes,
    n_z
):
    coolant_layers = [[] for _ in range(n_z)]

    for base_cell in coolant_cells:

        base_region = base_cell.region
        base_mat = base_cell.fill
        base_name = base_cell.name

        slice_regions = []

        # --- build axial slices
        for k in range(n_z):

            region = base_region & +z_planes[k] & -z_planes[k+1]
            slice_regions.append(region)

            c = openmc.Cell(
                fill=base_mat,
                region=region,
                name=f"{base_name}_z{k}"
            )

            universe.add_cell(c)
            coolant_layers[k].append(c)

        # --- union of slices
        sliced_union = slice_regions[0]
        for r in slice_regions[1:]:
            sliced_union |= r

        # --- remainder (top/bottom coolant)
        remainder_region = base_region & ~sliced_union

        remainder_cell = openmc.Cell(
            fill=base_mat,
            region=remainder_region,
            name=f"{base_name}_remainder"
        )

        # universe.add_cell(remainder_cell)

        # --- lower coolant (below heated region)

        lower_region = base_region & -z_planes[0]

        lower_cell = openmc.Cell(
            fill=base_mat,
            region=lower_region,
            name=f"{base_name}_lower"
        )
        universe.add_cell(lower_cell)

        # --- upper coolant (above heated region)
        upper_region = base_region & +z_planes[n_z]

        upper_cell = openmc.Cell(
            fill=base_mat,
            region=upper_region,
            name=f"{base_name}_upper"
        )
        universe.add_cell(upper_cell)

        # --- remove original coolant cell
        if base_cell.id in universe.cells:
            del universe.cells[base_cell.id]

    return coolant_layers, lower_cell, upper_cell




