# %%
import gmsh
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import root
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import glob
import re
import pandas as pd

# %% [markdown]
# # REGULAR 16-ROD FUEL ASSEMBLY (C7, D7, E7, F7, B6, B5, B4, F4, E3)

# %%


# ---------------------------------- meshing and stuff ---------------------------------------------------------------

# Initialize gmsh
gmsh.initialize()
gmsh.model.add('Regular_FE')

# Lets use opencascade cad kernel, to make the boolean operaitons easier
occ = gmsh.model.occ

r_matrix =  0.0035
r_al_clad = 0.0050


# The [0,0] coordinate is at the centre of FE
Fe_centres = [(-0.0255,-0.0255),
              (-0.0255,-0.0085),
              (-0.0255, 0.0085),
              (-0.0255, 0.0255),
              (-0.0085,-0.0255),
              (-0.0085,-0.0085),
              (-0.0085, 0.0085),
              (-0.0085, 0.0255),
              (0.0085,-0.0255),
              (0.0085,-0.0085),
              (0.0085, 0.0085),
              (0.0085, 0.0255),
              (0.0255,-0.0255),
              (0.0255,-0.0085),
              (0.0255, 0.0085),
              (0.0255, 0.0255)]


unit_cell_size = 0.072
box_channel_out = 0.068
box_channel_in = 0.065 


# --- 3. CREATE ALL OVERLAPPING PRIMITIVES ---
# We just throw all the shapes down on top of each other.
fluid_box = occ.addRectangle(-unit_cell_size/2, -unit_cell_size/2, 0, unit_cell_size, unit_cell_size)
clad_out = occ.addRectangle(-box_channel_out/2,-box_channel_out/2,0,box_channel_out,box_channel_out)
clad_in = occ.addRectangle(-box_channel_in/2,-box_channel_in/2,0,box_channel_in,box_channel_in)


clad_disks = []
fuel_disks = []
for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_al_clad, r_al_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_matrix, r_matrix))

# --- 4. THE MAGIC FRAGMENT OPERATION ---
# This shatters all overlapping shapes into distinct, non-overlapping pieces.
all_tools = [(2, tag) for tag in [clad_out,clad_in]+clad_disks + fuel_disks]
occ.fragment([(2, fluid_box)], all_tools)


# Must synchronize before we can ask gmsh for bounding boxes
occ.synchronize()

# --- 5. AUTOMATICALLY SORT THE SHATTERED PIECES ---
# Get all the new 2D surfaces created by the fragmentation
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box= None
fuel_tags = {i:[] for i in range(len(Fe_centres))}
clad_tags = {i:[] for i in range(len(Fe_centres))}




for dim, tag in surfaces:
    # get bounding boxes for each piece
    xmin, ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax-xmin,5)
    
    if width>=0.071:
        tag_out_water = tag
    elif width>=0.067:
        tag_al_box = tag
    elif width>=0.063:
        tag_in_water = tag
    else:
        com_x = 0.5*(xmin+xmax)
        com_y = 0.5*(ymin+ymax)
        
        for i, (cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)
            if dist<=0.001:
                if width>=2*r_al_clad:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)
    
    

gmsh.model.addPhysicalGroup(2,[tag_in_water], 300, name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water], 301, name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')


for i in range(len(Fe_centres)):
    pin_num = i + 1
    gmsh.model.addPhysicalGroup(2, fuel_tags[i], 100 + pin_num, name=f"Fuel_Pin_{pin_num}")
    gmsh.model.addPhysicalGroup(2, clad_tags[i], 200 + pin_num, name=f"Clad_Pin_{pin_num}")

# --- 7. GENERATE MESH ---
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0001)

print("Generating 2D Mesh...")
gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_4x4_fixed.msh")
print("Saved as 'bme_assembly_4x4_fixed.msh'")

# Visualize!
# gmsh.fltk.run()
# gmsh.finalize()



# %%
print("Number of nodes:", gmsh.model.mesh.getNodes()[0].size)
print("Number of elements:", len(gmsh.model.mesh.getElements()[1][0]))


# %% [markdown]
# #############################################################################################################################################

# %% [markdown]
# # REGULAR 13-ROD FUEL ASSEMBLY (B3)

# %%


# ---------------------------------- meshing and stuff ---------------------------------------------------------------

# Initialize gmsh
gmsh.initialize()
gmsh.model.add('B3_FE')

# Lets use opencascade cad kernel, to make the boolean operaitons easier
occ = gmsh.model.occ

r_matrix =  0.0035
r_al_clad = 0.0050


# The [0,0] coordinate is at the centre of FE
Fe_centres = [
              #(-0.0255,-0.0255),
              #(-0.0255,-0.0085),
              (-0.0255, 0.0085),
              (-0.0255, 0.0255),
              #(-0.0085,-0.0255),
              (-0.0085,-0.0085),
              (-0.0085, 0.0085),
              (-0.0085, 0.0255),
              (0.0085,-0.0255),
              (0.0085,-0.0085),
              (0.0085, 0.0085),
              (0.0085, 0.0255),
              (0.0255,-0.0255),
              (0.0255,-0.0085),
              (0.0255, 0.0085),
              (0.0255, 0.0255)
              ]

 

unit_cell_size = 0.072
box_channel_out = 0.068
box_channel_in = 0.065 


# --- 3. CREATE ALL OVERLAPPING PRIMITIVES ---
# We just throw all the shapes down on top of each other.
fluid_box = occ.addRectangle(-unit_cell_size/2, -unit_cell_size/2, 0, unit_cell_size, unit_cell_size)
clad_out = occ.addRectangle(-box_channel_out/2,-box_channel_out/2,0,box_channel_out,box_channel_out)
clad_in = occ.addRectangle(-box_channel_in/2,-box_channel_in/2,0,box_channel_in,box_channel_in)


clad_disks = []
fuel_disks = []
for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_al_clad, r_al_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_matrix, r_matrix))

# --- 4. THE MAGIC FRAGMENT OPERATION ---
# This shatters all overlapping shapes into distinct, non-overlapping pieces.
all_tools = [(2, tag) for tag in [clad_out,clad_in]+clad_disks + fuel_disks]
occ.fragment([(2, fluid_box)], all_tools)


# Must synchronize before we can ask gmsh for bounding boxes
occ.synchronize()

# --- 5. AUTOMATICALLY SORT THE SHATTERED PIECES ---
# Get all the new 2D surfaces created by the fragmentation
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box= None
fuel_tags = {i:[] for i in range(len(Fe_centres))}
clad_tags = {i:[] for i in range(len(Fe_centres))}




for dim, tag in surfaces:
    # get bounding boxes for each piece
    xmin, ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax-xmin,5)
    
    if width>=0.071:
        tag_out_water = tag
    elif width>=0.067:
        tag_al_box = tag
    elif width>=0.063:
        tag_in_water = tag
    else:
        com_x = 0.5*(xmin+xmax)
        com_y = 0.5*(ymin+ymax)
        
        for i, (cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)
            if dist<=0.001:
                if width>=0.011:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)
    
    

gmsh.model.addPhysicalGroup(2,[tag_in_water], 300, name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water], 301, name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')


for i in range(len(Fe_centres)):
    pin_num = i + 1
    gmsh.model.addPhysicalGroup(2, fuel_tags[i], 100 + pin_num, name=f"Fuel_Pin_{pin_num}")
    gmsh.model.addPhysicalGroup(2, clad_tags[i], 200 + pin_num, name=f"Clad_Pin_{pin_num}")

# --- 7. GENERATE MESH ---
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0001)

print("Generating 2D Mesh...")
gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_B3_fixed.msh")
print("Saved as 'bme_assembly_B3_fixed.msh'")

# Visualize!
# gmsh.fltk.run()
# gmsh.finalize()



# %% [markdown]
# ######################################################################################################################################################################################

# %% [markdown]
# # F3 FUEL ASSEMBLY 

# %%

import gmsh
import math

gmsh.initialize()
gmsh.model.add("assembly_12rod_irradiation_channel")

occ = gmsh.model.occ

# -----------------------------
# Geometry parameters [m]
# -----------------------------
r_fuel = 0.0035
r_clad = 0.0050

R_IR_INNER = 0.015
R_IR_OUTER = 0.016

assembly_pitch = 0.072
unit_cell_size = 0.072

box_channel_out = 0.068
box_channel_in  = 0.065

# -----------------------------
# Fuel rod centers [m]
# -----------------------------
Fe_centres = [
(-0.0255, -0.0255),
(-0.0255, -0.0085),
(-0.0255,  0.0085),
(-0.0255,  0.0255),
(-0.0085, -0.0255),
(-0.0085,  0.0255),
( 0.0085, -0.0255),
( 0.0085,  0.0255),
( 0.0255, -0.0255),
( 0.0255, -0.0085),
( 0.0255,  0.0085),
( 0.0255,  0.0255),
]

# -----------------------------
# Irradiation channel
# -----------------------------
ir_inner = occ.addDisk(0.0, 0.0, 0, R_IR_INNER, R_IR_INNER)
ir_outer = occ.addDisk(0.0, 0.0, 0, R_IR_OUTER, R_IR_OUTER)

# -----------------------------
# Assembly geometry
# -----------------------------
fluid_box = occ.addRectangle(-unit_cell_size/2, -unit_cell_size/2, 0,
                              unit_cell_size, unit_cell_size)

clad_out = occ.addRectangle(-box_channel_out/2, -box_channel_out/2, 0,
                             box_channel_out, box_channel_out)

clad_in = occ.addRectangle(-box_channel_in/2, -box_channel_in/2, 0,
                            box_channel_in, box_channel_in)

# -----------------------------
# Fuel rods
# -----------------------------
clad_disks = []
fuel_disks = []

for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_clad, r_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_fuel, r_fuel))

# -----------------------------
# Fragment everything
# -----------------------------
all_tools = [(2, tag) for tag in [clad_out, clad_in, ir_inner, ir_outer] + clad_disks + fuel_disks]

occ.fragment([(2, fluid_box)], all_tools)

occ.synchronize()

# -----------------------------
# Get resulting surfaces
# -----------------------------
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box = None
ir_inner_tag = None
ir_wall_tag = None

fuel_tags = {i: [] for i in range(len(Fe_centres))}
clad_tags = {i: [] for i in range(len(Fe_centres))}

tol = 1e-4

for dim, tag in surfaces:

    xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax - xmin,5)

    com_x = 0.5*(xmin+xmax)
    com_y = 0.5*(ymin+ymax)


    # irradiation channel
    dist_center = math.sqrt(com_x**2 + com_y**2)


    # outer coolant
    if width >= 0.071:
        tag_out_water = tag

    # aluminium cassette
    elif width >= 0.067:
        tag_al_box = tag

    # inner coolant
    elif width >= 0.063:
        tag_in_water = tag
    
    
    elif dist_center < R_IR_OUTER:
        if width >= 2*R_IR_OUTER :
            ir_wall_tag = tag
        else:
            ir_inner_tag = tag

    else:
        # rod regions
        for i,(cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)

            if dist<=0.001:
                if width>=2*r_clad:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)

# -----------------------------
# Physical groups
# -----------------------------
gmsh.model.addPhysicalGroup(2,[tag_in_water],300,name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water],301,name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')

gmsh.model.addPhysicalGroup(2,[ir_inner_tag],500,name='Irradiation_Water')
gmsh.model.addPhysicalGroup(2,[ir_wall_tag],501,name='Irradiation_Wall')

for i in range(len(Fe_centres)):
    pin = i+1
    gmsh.model.addPhysicalGroup(2,fuel_tags[i],100+pin,name=f"Fuel_Pin_{pin}")
    gmsh.model.addPhysicalGroup(2,clad_tags[i],200+pin,name=f"Clad_Pin_{pin}")

# -----------------------------
# Mesh
# -----------------------------
gmsh.option.setNumber("Mesh.MeshSizeMax",0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin",0.0001)

gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_F3_fixed.msh")

# gmsh.fltk.run()
# gmsh.finalize()

 

# %% [markdown]
# #############################################################################################################################################################################################

# %% [markdown]
# # LOWER-RIGHT CUT (F5, E4, D3)

# %%

import gmsh
import math

gmsh.initialize()
gmsh.model.add("assembly_lower_right_cut")

occ = gmsh.model.occ

# -----------------------------
# Geometry parameters [m]
# -----------------------------
r_fuel = 0.0035
r_clad = 0.0050

unit_cell_size = 0.072

box_channel_out = 0.068
box_channel_in  = 0.065
# -----------------------------
# Fuel rod centers [m]
# -----------------------------

# x0 = 0.0
# y0 = 10.8

# Fe_centres_shifted = [(x - x0, y - y0) for x, y in Fe_centres]

Fe_centres = [
(-0.0255, -0.0255),
(-0.0255, -0.0085),
(-0.0255,  0.0085),
(-0.0255,  0.0255),

(-0.0105, -0.0255),
(-0.0085, -0.0085),
(-0.0085,  0.0085),
(-0.0085,  0.0255),

(0.0045, -0.0255),
(0.0085, -0.0085),
(0.0085,  0.0085),
(0.0085,  0.0255),

(0.0180, -0.0180),
(0.0255, -0.0045),
(0.0255,  0.0105),
(0.0255,  0.0255)
]

# Fe_centres = [(x/100, y/100) for x, y in Fe_centres]

box_outer_points = [
(-0.034, -0.034),
(-0.034,  0.034),
( 0.034,  0.034),
( 0.034, -0.014),
( 0.014, -0.034),
(-0.034, -0.034)
]


box_inner_points = [
(-0.0325, -0.0325),
(-0.0325,  0.0325),
( 0.0325,  0.0325),
( 0.0325, -0.0135),
( 0.0135, -0.0325),
(-0.0325, -0.0325)
]

p = []
for x, y in box_outer_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_out = occ.addPlaneSurface([loop])

p = []
for x, y in box_inner_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_in = occ.addPlaneSurface([loop])


half = 0.036
r_cut = 0.0145
xc, yc = 0.036, -0.036

p1 = occ.addPoint(-0.036, -0.036, 0)
p2 = occ.addPoint(-0.036,  0.036, 0)
p3 = occ.addPoint( 0.036,  0.036, 0)
p4 = occ.addPoint( 0.036, -0.0215, 0)
p5 = occ.addPoint( 0.0215, -0.036, 0)
pc = occ.addPoint( 0.036, -0.036, 0)


# lines
l1 = occ.addLine(p1, p2)
l2 = occ.addLine(p2, p3)
l3 = occ.addLine(p3, p4)
l4 = occ.addCircleArc(p4, pc, p5)
l5 = occ.addLine(p5, p1)

loop = occ.addCurveLoop([l1, l2, l3, l4, l5])
fluid_box = occ.addPlaneSurface([loop])

# -----------------------------
# Fuel rods
# -----------------------------
clad_disks = []
fuel_disks = []

for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_clad, r_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_fuel, r_fuel))

# -----------------------------
# Fragment everything
# -----------------------------
all_tools = [(2, tag) for tag in [box_out, box_in] + clad_disks + fuel_disks]

occ.fragment([(2, fluid_box)], all_tools)

occ.synchronize()

# -----------------------------
# Get resulting surfaces
# -----------------------------
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box = None


fuel_tags = {i: [] for i in range(len(Fe_centres))}
clad_tags = {i: [] for i in range(len(Fe_centres))}

tol = 1e-4

for dim, tag in surfaces:

    xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax - xmin,5)

    com_x = 0.5*(xmin+xmax)
    com_y = 0.5*(ymin+ymax)


    # irradiation channel
    dist_center = math.sqrt(com_x**2 + com_y**2)


    # outer coolant
    if width >= 0.071:
        tag_out_water = tag

    # aluminium cassette
    elif width >= 0.067:
        tag_al_box = tag

    # inner coolant
    elif width >= 0.063:
        tag_in_water = tag


    else:
        # rod regions
        for i,(cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)

            if dist<=0.001:
                if width>=2*r_clad:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)

# -----------------------------
# Physical groups
# -----------------------------
gmsh.model.addPhysicalGroup(2,[tag_in_water],300,name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water],301,name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')

for i in range(len(Fe_centres)):
    pin = i+1
    gmsh.model.addPhysicalGroup(2,fuel_tags[i],100+pin,name=f"Fuel_Pin_{pin}")
    gmsh.model.addPhysicalGroup(2,clad_tags[i],200+pin,name=f"Clad_Pin_{pin}")

# -----------------------------
# Mesh
# -----------------------------
gmsh.option.setNumber("Mesh.MeshSizeMax",0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin",0.0001)

gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_F5_fixed.msh")

# gmsh.fltk.run()
# gmsh.finalize()

 

# %% [markdown]
# ######################################################################################################################################################################################

# %% [markdown]
# # LOWER-LEFT CUT (F6, D6)

# %%

import gmsh
import math

gmsh.initialize()
gmsh.model.add("assembly_lower_left")

occ = gmsh.model.occ

# -----------------------------
# Geometry parameters [m]
# -----------------------------
r_fuel = 0.0035
r_clad = 0.0050

unit_cell_size = 0.072

box_channel_out = 0.068
box_channel_in  = 0.065
# -----------------------------
# Fuel rod centers [m]
# -----------------------------

# x0 = 0.0
# y0 = 10.8

# Fe_centres_shifted = [(x - x0, y - y0) for x, y in Fe_centres]

Fe_centres = [
(-0.0255, -0.0255),
(-0.0255, -0.0085),
(-0.0255,  0.0085),
(-0.0255,  0.0255),

(-0.0105, -0.0255),
(-0.0085, -0.0085),
(-0.0085,  0.0085),
(-0.0085,  0.0255),

(0.0045, -0.0255),
(0.0085, -0.0085),
(0.0085,  0.0085),
(0.0085,  0.0255),

(0.0180, -0.0180),
(0.0255, -0.0045),
(0.0255,  0.0105),
(0.0255,  0.0255)
]


# Fe_centres = [(x/100, y/100) for x, y in Fe_centres]

box_outer_points = [
(-0.034, -0.034),
(-0.034,  0.034),
( 0.034,  0.034),
( 0.034, -0.014),
( 0.014, -0.034),
(-0.034, -0.034)
]



box_inner_points = [
(-0.0325, -0.0325),
(-0.0325,  0.0325),
( 0.0325,  0.0325),
( 0.0325, -0.0135),
( 0.0135, -0.0325),
(-0.0325, -0.0325)
]


def mirror_points(points):
    return [(-x, y) for (x, y) in points]


Fe_centres = mirror_points(Fe_centres)
box_outer_points = mirror_points(box_outer_points)
box_inner_points = mirror_points(box_inner_points)

p = []
for x, y in box_outer_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_out = occ.addPlaneSurface([loop])

p = []
for x, y in box_inner_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_in = occ.addPlaneSurface([loop])



half = 0.036
r_cut = 0.0145
xc, yc =- 0.036, -0.036

p1 = occ.addPoint( 0.036, -0.036, 0)
p2 = occ.addPoint( 0.036,  0.036, 0)
p3 = occ.addPoint(-0.036,  0.036, 0)
p4 = occ.addPoint(-0.036, -0.0215, 0)
p5 = occ.addPoint(-0.0215, -0.036, 0)
pc = occ.addPoint(-0.036, -0.036, 0)


# lines
l1 = occ.addLine(p1, p2)
l2 = occ.addLine(p2, p3)
l3 = occ.addLine(p3, p4)
l4 = occ.addCircleArc(p4, pc, p5)
l5 = occ.addLine(p5, p1)

loop = occ.addCurveLoop([l1, l2, l3, l4, l5])
fluid_box = occ.addPlaneSurface([loop])

# -----------------------------
# Fuel rods
# -----------------------------
clad_disks = []
fuel_disks = []

for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_clad, r_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_fuel, r_fuel))

# -----------------------------
# Fragment everything
# -----------------------------
all_tools = [(2, tag) for tag in [box_out, box_in] + clad_disks + fuel_disks]

occ.fragment([(2, fluid_box)], all_tools)

occ.synchronize()

# -----------------------------
# Get resulting surfaces
# -----------------------------
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box = None


fuel_tags = {i: [] for i in range(len(Fe_centres))}
clad_tags = {i: [] for i in range(len(Fe_centres))}

tol = 1e-4

for dim, tag in surfaces:

    xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax - xmin,5)

    com_x = 0.5*(xmin+xmax)
    com_y = 0.5*(ymin+ymax)


    # irradiation channel
    dist_center = math.sqrt(com_x**2 + com_y**2)


    # outer coolant
    if width >= 0.071:
        tag_out_water = tag

    # aluminium cassette
    elif width >= 0.067:
        tag_al_box = tag

    # inner coolant
    elif width >= 0.063:
        tag_in_water = tag


    else:
        # rod regions
        for i,(cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)

            if dist<=0.001:
                if width>=2*r_clad:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)

# -----------------------------
# Physical groups
# -----------------------------
gmsh.model.addPhysicalGroup(2,[tag_in_water],300,name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water],301,name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')

for i in range(len(Fe_centres)):
    pin = i+1
    gmsh.model.addPhysicalGroup(2,fuel_tags[i],100+pin,name=f"Fuel_Pin_{pin}")
    gmsh.model.addPhysicalGroup(2,clad_tags[i],200+pin,name=f"Clad_Pin_{pin}")

# -----------------------------
# Mesh
# -----------------------------
gmsh.option.setNumber("Mesh.MeshSizeMax",0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin",0.0001)

gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_F6_fixed.msh")

# gmsh.fltk.run()
# gmsh.finalize()

 

# %% [markdown]
# ##############################################################################################################################################################################################

# %% [markdown]
# # UPPER-LEFT CUT (E6, C4, C6)

# %%

import gmsh
import math

gmsh.initialize()
gmsh.model.add("assembly_upper_left")

occ = gmsh.model.occ

# -----------------------------
# Geometry parameters [m]
# -----------------------------
r_fuel = 0.0035
r_clad = 0.0050

unit_cell_size = 0.072

box_channel_out = 0.068
box_channel_in  = 0.065
# -----------------------------
# Fuel rod centers [m]
# -----------------------------

# x0 = 0.0
# y0 = 10.8

# Fe_centres_shifted = [(x - x0, y - y0) for x, y in Fe_centres]

Fe_centres = [
(-0.0255, -0.0255),
(-0.0255, -0.0085),
(-0.0255,  0.0085),
(-0.0255,  0.0255),

(-0.0105, -0.0255),
(-0.0085, -0.0085),
(-0.0085,  0.0085),
(-0.0085,  0.0255),

(0.0045, -0.0255),
(0.0085, -0.0085),
(0.0085,  0.0085),
(0.0085,  0.0255),

(0.0180, -0.0180),
(0.0255, -0.0045),
(0.0255,  0.0105),
(0.0255,  0.0255)
]


# Fe_centres = [(x/100, y/100) for x, y in Fe_centres]

box_outer_points = [
(-0.034, -0.034),
(-0.034,  0.034),
( 0.034,  0.034),
( 0.034, -0.014),
( 0.014, -0.034),
(-0.034, -0.034)
]



box_inner_points = [
(-0.0325, -0.0325),
(-0.0325,  0.0325),
( 0.0325,  0.0325),
( 0.0325, -0.0135),
( 0.0135, -0.0325),
(-0.0325, -0.0325)
]


def mirror_points(points):
    return [(-x, -y) for (x, y) in points]


Fe_centres = mirror_points(Fe_centres)
box_outer_points = mirror_points(box_outer_points)
box_inner_points = mirror_points(box_inner_points)

p = []
for x, y in box_outer_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_out = occ.addPlaneSurface([loop])

p = []
for x, y in box_inner_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_in = occ.addPlaneSurface([loop])



half = 0.036
r_cut = 0.0145
xc, yc = -0.036, 0.036



p1 = occ.addPoint( 0.036,  0.036, 0)
p2 = occ.addPoint( 0.036, -0.036, 0)
p3 = occ.addPoint(-0.036, -0.036, 0)
p4 = occ.addPoint(-0.036,  0.0215, 0)
p5 = occ.addPoint(-0.0215,  0.036, 0)
pc = occ.addPoint(-0.036,  0.036, 0)


# lines
l1 = occ.addLine(p1, p2)
l2 = occ.addLine(p2, p3)
l3 = occ.addLine(p3, p4)
l4 = occ.addCircleArc(p5, pc, p4)
l5 = occ.addLine(p5, p1)

loop = occ.addCurveLoop([l1, l2, l3, l4, l5])
fluid_box = occ.addPlaneSurface([loop])

# -----------------------------
# Fuel rods
# -----------------------------
clad_disks = []
fuel_disks = []

for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_clad, r_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_fuel, r_fuel))

# -----------------------------
# Fragment everything
# -----------------------------
all_tools = [(2, tag) for tag in [box_out, box_in] + clad_disks + fuel_disks]

occ.fragment([(2, fluid_box)], all_tools)

occ.synchronize()

# -----------------------------
# Get resulting surfaces
# -----------------------------
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box = None


fuel_tags = {i: [] for i in range(len(Fe_centres))}
clad_tags = {i: [] for i in range(len(Fe_centres))}

tol = 1e-4

for dim, tag in surfaces:

    xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax - xmin,5)

    com_x = 0.5*(xmin+xmax)
    com_y = 0.5*(ymin+ymax)


    # irradiation channel
    dist_center = math.sqrt(com_x**2 + com_y**2)


    # outer coolant
    if width >= 0.071:
        tag_out_water = tag

    # aluminium cassette
    elif width >= 0.067:
        tag_al_box = tag

    # inner coolant
    elif width >= 0.063:
        tag_in_water = tag


    else:
        # rod regions
        for i,(cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)

            if dist<=0.001:
                if width>=2*r_clad:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)

# -----------------------------
# Physical groups
# -----------------------------
gmsh.model.addPhysicalGroup(2,[tag_in_water],300,name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water],301,name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')

for i in range(len(Fe_centres)):
    pin = i+1
    gmsh.model.addPhysicalGroup(2,fuel_tags[i],100+pin,name=f"Fuel_Pin_{pin}")
    gmsh.model.addPhysicalGroup(2,clad_tags[i],200+pin,name=f"Clad_Pin_{pin}")

# -----------------------------
# Mesh
# -----------------------------
gmsh.option.setNumber("Mesh.MeshSizeMax",0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin",0.0001)

gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_E6_fixed.msh")

# gmsh.fltk.run()
# gmsh.finalize()

 

# %% [markdown]
# #############################################################################################################################################################################################

# %% [markdown]
# # UPPER-RIGHT CUT (C3, C5)

# %%

import gmsh
import math

gmsh.initialize()
gmsh.model.add("assembly_upper_right")

occ = gmsh.model.occ

# -----------------------------
# Geometry parameters [m]
# -----------------------------
r_fuel = 0.0035
r_clad = 0.0050

unit_cell_size = 0.072

box_channel_out = 0.068
box_channel_in  = 0.065
# -----------------------------
# Fuel rod centers [m]
# -----------------------------

# x0 = 0.0
# y0 = 10.8

# Fe_centres_shifted = [(x - x0, y - y0) for x, y in Fe_centres]

Fe_centres = [
(-0.0255, -0.0255),
(-0.0255, -0.0085),
(-0.0255,  0.0085),
(-0.0255,  0.0255),

(-0.0105, -0.0255),
(-0.0085, -0.0085),
(-0.0085,  0.0085),
(-0.0085,  0.0255),

(0.0045, -0.0255),
(0.0085, -0.0085),
(0.0085,  0.0085),
(0.0085,  0.0255),

(0.0180, -0.0180),
(0.0255, -0.0045),
(0.0255,  0.0105),
(0.0255,  0.0255)
]


# Fe_centres = [(x/100, y/100) for x, y in Fe_centres]

box_outer_points = [
(-0.034, -0.034),
(-0.034,  0.034),
( 0.034,  0.034),
( 0.034, -0.014),
( 0.014, -0.034),
(-0.034, -0.034)
]



box_inner_points = [
(-0.0325, -0.0325),
(-0.0325,  0.0325),
( 0.0325,  0.0325),
( 0.0325, -0.0135),
( 0.0135, -0.0325),
(-0.0325, -0.0325)
]


def mirror_points(points):
    return [(x, -y) for (x, y) in points]


Fe_centres = mirror_points(Fe_centres)
box_outer_points = mirror_points(box_outer_points)
box_inner_points = mirror_points(box_inner_points)

p = []
for x, y in box_outer_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_out = occ.addPlaneSurface([loop])

p = []
for x, y in box_inner_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_in = occ.addPlaneSurface([loop])



half = 0.036
r_cut = 0.0145
xc, yc = 0.036, 0.036

p1 = occ.addPoint(-half, -half, 0)
p2 = occ.addPoint( half, -half, 0)
p3 = occ.addPoint( half,  half-r_cut, 0)
p4 = occ.addPoint( half-r_cut,  half, 0)
p5 = occ.addPoint(-half,  half, 0)
pc = occ.addPoint( half,  half, 0)

# lines
l1 = occ.addLine(p1, p2)
l2 = occ.addLine(p2, p3)
l3 = occ.addCircleArc(p3, pc, p4)
l4 = occ.addLine(p4, p5)
l5 = occ.addLine(p5, p1)

loop = occ.addCurveLoop([l1, l2, l3, l4, l5])
fluid_box = occ.addPlaneSurface([loop])

# -----------------------------
# Fuel rods
# -----------------------------
clad_disks = []
fuel_disks = []

for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_clad, r_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_fuel, r_fuel))

# -----------------------------
# Fragment everything
# -----------------------------
all_tools = [(2, tag) for tag in [box_out, box_in] + clad_disks + fuel_disks]

occ.fragment([(2, fluid_box)], all_tools)

occ.synchronize()

# -----------------------------
# Get resulting surfaces
# -----------------------------
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box = None


fuel_tags = {i: [] for i in range(len(Fe_centres))}
clad_tags = {i: [] for i in range(len(Fe_centres))}

tol = 1e-4

for dim, tag in surfaces:

    xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax - xmin,5)

    com_x = 0.5*(xmin+xmax)
    com_y = 0.5*(ymin+ymax)


    # irradiation channel
    dist_center = math.sqrt(com_x**2 + com_y**2)


    # outer coolant
    if width >= 0.071:
        tag_out_water = tag

    # aluminium cassette
    elif width >= 0.067:
        tag_al_box = tag

    # inner coolant
    elif width >= 0.063:
        tag_in_water = tag


    else:
        # rod regions
        for i,(cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)

            if dist<=0.001:
                if width>=2*r_clad:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)

# -----------------------------
# Physical groups
# -----------------------------
gmsh.model.addPhysicalGroup(2,[tag_in_water],300,name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water],301,name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')

for i in range(len(Fe_centres)):
    pin = i+1
    gmsh.model.addPhysicalGroup(2,fuel_tags[i],100+pin,name=f"Fuel_Pin_{pin}")
    gmsh.model.addPhysicalGroup(2,clad_tags[i],200+pin,name=f"Clad_Pin_{pin}")

# -----------------------------
# Mesh
# -----------------------------
gmsh.option.setNumber("Mesh.MeshSizeMax",0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin",0.0001)

gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_C3_fixed.msh")

# gmsh.fltk.run()
# gmsh.finalize()

 

# %% [markdown]
# #############################################################################################################################################################################################

# %% [markdown]
# # 3-CORNER-CUT ASSEMBLIES (D4, E5)

# %%

import gmsh
import math

gmsh.initialize()
gmsh.model.add("assembly_three_corner_cut")

occ = gmsh.model.occ

# -----------------------------
# Geometry parameters [m]
# -----------------------------
r_fuel = 0.0035
r_clad = 0.0050

unit_cell_size = 0.072

box_channel_out = 0.068
box_channel_in  = 0.065
# -----------------------------
# Fuel rod centers [m]
# -----------------------------

# x0 = 0.0
# y0 = 10.8

# Fe_centres_shifted = [(x - x0, y - y0) for x, y in Fe_centres]

Fe_centres = [

(0.0255,  0.0035),
(0.0255, -0.0110),
(0.0255, -0.0255),
(0.0190,  0.0165),

(0.0110, -0.0255),
(0.0085,  0.0085),
(0.0085, -0.0085),
(0.0070,  0.0255),

(-0.0035, -0.0255),
(-0.0075,  0.0255),
(-0.0085,  0.0085),
(-0.0085, -0.0085),

(-0.0185,  0.0185),
(-0.0165, -0.0190),
(-0.0255,  0.0075),
(-0.0255, -0.0070)

]
# Fe_centres = [(x/100, y/100) for x, y in Fe_centres]

box_outer_points = [
    ( 0.034, -0.034),
    (-0.014, -0.034),
    (-0.034, -0.014),
    (-0.034,  0.014),
    (-0.014,  0.034),
    ( 0.014,  0.034),
    ( 0.034,  0.014),
    ( 0.034, -0.034)
]

box_inner_points = [
( 0.0325, -0.0325),
(-0.0135, -0.0325),
(-0.0325, -0.0135),
(-0.0325,  0.0135),
(-0.0135,  0.0325),
( 0.0135,  0.0325),
( 0.0325,  0.0135),
( 0.0325, -0.0325)
]



p = []
for x, y in box_outer_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_out = occ.addPlaneSurface([loop])

p = []
for x, y in box_inner_points:
    p.append(occ.addPoint(x, y, 0))

lines = []
for i in range(len(p)-1):
    lines.append(occ.addLine(p[i], p[i+1]))

loop = occ.addCurveLoop(lines)
box_in = occ.addPlaneSurface([loop])


def make_two_arc_box(half, r_cut):
    """
    Creates a surface of a square-like box with:
      - lower-left corner replaced by a quarter circle
      - upper-right corner replaced by a quarter circle
    Returns surface tag.
    """

    # centers of the two quarter-circles
    c_ll = occ.addPoint(-half, -half, 0)
    c_ur = occ.addPoint( half,  half, 0)

    # boundary points, counter-clockwise
    p1 = occ.addPoint(-half,      -half + r_cut, 0)   # left side, start LL arc
    p2 = occ.addPoint(-half + r_cut, -half,      0)   # bottom side, end LL arc
    p3 = occ.addPoint( half,      -half,         0)   # bottom-right
    p4 = occ.addPoint( half,       half - r_cut, 0)   # right side, start UR arc
    p5 = occ.addPoint( half - r_cut, half,       0)   # top side, end UR arc
    p6 = occ.addPoint(-half,       half,         0)   # top-left

    # curves
    l1 = occ.addCircleArc(p1, c_ll, p2)   # lower-left quarter circle
    l2 = occ.addLine(p2, p3)              # bottom edge
    l3 = occ.addLine(p3, p4)              # right edge
    l4 = occ.addCircleArc(p4, c_ur, p5)   # upper-right quarter circle
    l5 = occ.addLine(p5, p6)              # top edge
    l6 = occ.addLine(p6, p1)              # left edge

    loop = occ.addCurveLoop([l1, l2, l3, l4, l5, l6])
    return occ.addPlaneSurface([loop])

half = 0.036
r_cut = 0.0145
xc, yc = 0.036, 0.036

half = 0.036
r_cut = 0.0145

fluid_box = make_two_arc_box(half, r_cut)


# -----------------------------
# Fuel rods
# -----------------------------
clad_disks = []
fuel_disks = []

for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_clad, r_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_fuel, r_fuel))

# -----------------------------
# Fragment everything
# -----------------------------
all_tools = [(2, tag) for tag in [box_out, box_in] + clad_disks + fuel_disks]

occ.fragment([(2, fluid_box)], all_tools)

occ.synchronize()

# -----------------------------
# Get resulting surfaces
# -----------------------------
surfaces = occ.getEntities(2)

tag_in_water = None
tag_out_water = None
tag_al_box = None


fuel_tags = {i: [] for i in range(len(Fe_centres))}
clad_tags = {i: [] for i in range(len(Fe_centres))}

tol = 1e-4

for dim, tag in surfaces:

    xmin,ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
    width = round(xmax - xmin,5)

    com_x = 0.5*(xmin+xmax)
    com_y = 0.5*(ymin+ymax)


    # irradiation channel
    dist_center = math.sqrt(com_x**2 + com_y**2)


    # outer coolant
    if width >= 0.071:
        tag_out_water = tag

    # aluminium cassette
    elif width >= 0.067:
        tag_al_box = tag

    # inner coolant
    elif width >= 0.063:
        tag_in_water = tag


    else:
        # rod regions
        for i,(cx,cy) in enumerate(Fe_centres):
            dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)

            if dist<=0.001:
                if width>=2*r_clad:
                    clad_tags[i].append(tag)
                else:
                    fuel_tags[i].append(tag)

# -----------------------------
# Physical groups
# -----------------------------
gmsh.model.addPhysicalGroup(2,[tag_in_water],300,name='Inner_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_out_water],301,name='Outer_Coolant')
gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')

for i in range(len(Fe_centres)):
    pin = i+1
    gmsh.model.addPhysicalGroup(2,fuel_tags[i],100+pin,name=f"Fuel_Pin_{pin}")
    gmsh.model.addPhysicalGroup(2,clad_tags[i],200+pin,name=f"Clad_Pin_{pin}")

# -----------------------------
# Mesh
# -----------------------------
gmsh.option.setNumber("Mesh.MeshSizeMax",0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMin",0.0001)

gmsh.model.mesh.generate(2)
gmsh.write("bme_assembly_E5_fixed.msh")

# gmsh.fltk.run()
# gmsh.finalize()

 

# %% [markdown]
# #############################################################################################################################################################################################

# %% [markdown]
# # D5

# %%
import gmsh
import math

gmsh.initialize()
gmsh.model.add("bme_centered_air_chamber")
occ = gmsh.model.occ

# -------------------------------------------------
# parameters
# -------------------------------------------------
r_fuel = 0.0035
r_clad = 0.0050
r_cut = 0.0145  # outer fluid box corner cut radius [m]

# -------------------------------------------------
# helpers
# -------------------------------------------------
def cm(x, y, y_shift_cm=0.0):
    return (x * 1e-2, (y + y_shift_cm) * 1e-2)

# def add_polygon_surface(points):
#     pts = [occ.addPoint(x, y, 0) for x, y in points]
#     lines = [occ.addLine(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
#     loop = occ.addCurveLoop(lines)
#     return occ.addPlaneSurface([loop])

def point_in_polygon(x, y, poly):
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if ((y1 > y) != (y2 > y)):
            xinters = (x2 - x1) * (y - y1) / (y2 - y1 + 1e-30) + x1
            if x < xinters:
                inside = not inside
    return inside

# -------------------------------------------------
# cassette geometry
# -------------------------------------------------
def add_polygon_surface(points):
    pts = [occ.addPoint(x, y, 0) for x, y in points]
    lines = [occ.addLine(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    loop = occ.addCurveLoop(lines)
    return occ.addPlaneSurface([loop])

# -------------------------------------------------
# cassette geometry
# -------------------------------------------------
box_outer_points = [
    (-0.0340,  0.0140),
    (-0.0140,  0.0340),
    ( 0.0340,  0.0340),
    ( 0.0340, -0.0140),
    ( 0.0140, -0.0340),
    (-0.0140, -0.0340),
    (-0.0340, -0.0140),
    (-0.0340,  0.0140)
]

box_inner_points = [
    (-0.0325,  0.0135),
    (-0.0135,  0.0325),
    ( 0.0325,  0.0325),
    ( 0.0325, -0.0135),
    ( 0.0135, -0.0325),
    (-0.0135, -0.0325),
    (-0.0325, -0.0135),
    (-0.0325,  0.0135)
]

box_out = add_polygon_surface(box_outer_points)
box_in  = add_polygon_surface(box_inner_points)



# -------------------------------------------------
# air region from your cm points
# -------------------------------------------------


# -------------------------------------------------
# air chamber geometry
# -------------------------------------------------

p1 = occ.addPoint( 0.0070, -0.0308, 0)
p2 = occ.addPoint( 0.0070, -0.0213, 0)

p3 = occ.addPoint(-0.0005, -0.0080, 0)
p4 = occ.addPoint(-0.0005,  0.0080, 0)

p5 = occ.addPoint( 0.0070,  0.0213, 0)
p6 = occ.addPoint( 0.0070,  0.0308, 0)

p7 = occ.addPoint( 0.0308,  0.0308, 0)
p8 = occ.addPoint( 0.0308, -0.0130, 0)
p9 = occ.addPoint( 0.0130, -0.0308, 0)

# arc centers
c_lower = occ.addPoint(0.015, -0.008, 0)
c_upper = occ.addPoint(0.015,  0.008, 0)

# curves
l1 = occ.addLine(p1, p2)

arc1 = occ.addSpline([p2, p3])


l2 = occ.addLine(p3, p4)

arc2 = occ.addSpline([p4, p5])


l3 = occ.addLine(p5, p6)
l4 = occ.addLine(p6, p7)
l5 = occ.addLine(p7, p8)
l6 = occ.addLine(p8, p9)
l7 = occ.addLine(p9, p1)

air_loop = occ.addCurveLoop([l1, arc1, l2, arc2, l3, l4, l5, l6, l7])
air_region = occ.addPlaneSurface([air_loop])
# -------------------------------------------------
# separator wall from your OUTER wall points in cm
# -------------------------------------------------




p11 = occ.addPoint( 0.0050, -0.0325, 0)
p12 = occ.addPoint( 0.0050, -0.02236, 0)

p13 = occ.addPoint(-0.0025, -0.0080, 0)
p14 = occ.addPoint(-0.0025,  0.0080, 0)

p15 = occ.addPoint( 0.0050,  0.02236, 0)
p16 = occ.addPoint( 0.0050,  0.0325, 0)

p17 = occ.addPoint( 0.0325,  0.0325, 0)
p18 = occ.addPoint( 0.0325, -0.0135, 0)
p19 = occ.addPoint( 0.0135, -0.0325, 0)

# arc centers
c_lower = occ.addPoint(0.015, -0.008, 0)
c_upper = occ.addPoint(0.015,  0.008, 0)

# curves
l11 = occ.addLine(p11, p12)

arc11 = occ.addSpline([p12, p13])


l12 = occ.addLine(p13, p14)

arc12 = occ.addSpline([p14, p15])


l13 = occ.addLine(p15, p16)
l14 = occ.addLine(p16, p17)
l15 = occ.addLine(p17, p18)
l16 = occ.addLine(p18, p19)
l17 = occ.addLine(p19, p11)

air_wall = occ.addCurveLoop([l11, arc11, l12, arc12, l13, l14, l15, l16, l17])
air_wall_region = occ.addPlaneSurface([air_wall])


# -------------------------------------------------
# fuel rods
# -------------------------------------------------
Fe_centres = [
    (-0.0035,  0.0255),
    (-0.0075, -0.0255),
    (-0.0085,  0.0085),
    (-0.0085, -0.0085),
    (-0.0165,  0.0190),
    (-0.0185, -0.0185),
    (-0.0255,  0.0070),
    (-0.0255, -0.0075),
]

clad_disks = []
fuel_disks = []
for cx, cy in Fe_centres:
    clad_disks.append(occ.addDisk(cx, cy, 0, r_clad, r_clad))
    fuel_disks.append(occ.addDisk(cx, cy, 0, r_fuel, r_fuel))

# -------------------------------------------------
# outer coolant region (unit cell 0.072 m)
# upper-left and lower-right quarter-circle cuts
# -------------------------------------------------
r_cut = 0.0145

# square corners
xmin = -0.036
xmax =  0.036
ymin = -0.036
ymax =  0.036

# arc points
p_ul1 = occ.addPoint(xmin + r_cut, ymax, 0)
p_ul2 = occ.addPoint(xmin, ymax - r_cut, 0)

p_lr1 = occ.addPoint(xmax - r_cut, ymin, 0)
p_lr2 = occ.addPoint(xmax, ymin + r_cut, 0)

# remaining square points
p_tr = occ.addPoint(xmax, ymax, 0)
p_bl = occ.addPoint(xmin, ymin, 0)

# arc centers
c_ul = occ.addPoint(xmin, ymax, 0)
c_lr = occ.addPoint(xmax, ymin, 0)

# boundary curves
l1 = occ.addLine(p_ul1, p_tr)
l2 = occ.addLine(p_tr, p_lr2)

arc_lr = occ.addCircleArc(p_lr2, c_lr, p_lr1)

l3 = occ.addLine(p_lr1, p_bl)
l4 = occ.addLine(p_bl, p_ul2)

arc_ul = occ.addCircleArc(p_ul2, c_ul, p_ul1)

# loop
fluid_loop = occ.addCurveLoop([l1, l2, arc_lr, l3, l4, arc_ul])
fluid_box = occ.addPlaneSurface([fluid_loop])

# subtract cassette exterior
outer_fluid, _ = occ.cut(
    [(2, fluid_box)],
    [(2, box_out)],
    removeTool=False
)
# -------------------------------------------------
# fragment geometry
# -------------------------------------------------
objects = [(2, box_in), (2, box_out)] + outer_fluid
tools = (
    [(2, air_region), (2, air_wall_region)] +
    [(2, s) for s in clad_disks] +
    [(2, s) for s in fuel_disks]
)

occ.fragment(objects, tools)
occ.synchronize()

# -------------------------------------------------
# classify regions
# -------------------------------------------------
surfaces = gmsh.model.getEntities(2)

fuel_tags = {i: [] for i in range(len(Fe_centres))}
clad_tags = {i: [] for i in range(len(Fe_centres))}

tag_air = []
tag_wall = []
tag_coolant_inside = []
tag_coolant_outside = []
tag_can = []

# -------------------------------------------------
# classify regions
# -------------------------------------------------
# helper: test if any bbox corner is inside polygon
def bbox_any_corner_in_polygon(xmin, ymin, xmax, ymax, poly):
    return (
        point_in_polygon(xmin, ymin, poly) or
        point_in_polygon(xmax, ymin, poly) or
        point_in_polygon(xmax, ymax, poly) or
        point_in_polygon(xmin, ymax, poly)
    )


# classify regions
# -------------------------------------------------
for dim, tag in surfaces:
    xmin_b, ymin_b, zmin_b, xmax_b, ymax_b, zmax_b = gmsh.model.getBoundingBox(dim, tag)

    cx = 0.5 * (xmin_b + xmax_b)
    cy = 0.5 * (ymin_b + ymax_b)
    width  = xmax_b - xmin_b
    height = ymax_b - ymin_b

    # ------------------------------
    # rods
    # ------------------------------
    assigned = False
    for i, (fx, fy) in enumerate(Fe_centres):
        dist = math.sqrt((cx - fx)**2 + (cy - fy)**2)

        if dist <= 0.001:
            if width >= 2 * r_clad - 5e-4:
                clad_tags[i].append(tag)
            else:
                fuel_tags[i].append(tag)
            assigned = True
            break

    if assigned:
        continue

    # ------------------------------
    # large regions by width
    # ------------------------------
    w_cm = width * 100.0

    # outer coolant: full unit cell ~ 7.2 cm
    if w_cm > 7.0:
        tag_coolant_outside.append(tag)
        continue

    # assembly duct: ~ 6.8 cm
    if 6.6 < w_cm < 6.95:
        tag_can.append(tag)
        continue

    # inner coolant: left side region ~ 3.75 cm
    if 3.6 < w_cm < 3.9:
        tag_coolant_inside.append(tag)
        continue

    # air chamber: ~ 3.13 cm
    if 3.0 < w_cm < 3.3:
        tag_air.append(tag)
        continue

    # separator wall: narrow band around 2 mm, but fragmented pieces may be larger
    if xmin_b > -0.003 and xmax_b < 0.033:
        tag_wall.append(tag)
        continue

    # fallback
    tag_coolant_inside.append(tag)
# -------------------------------------------------
# physical groups
# -------------------------------------------------
if tag_coolant_inside:
    gmsh.model.addPhysicalGroup(2, tag_coolant_inside, 300)
    gmsh.model.setPhysicalName(2, 300, "Coolant_Inside")

if tag_air:
    gmsh.model.addPhysicalGroup(2, tag_air, 301)
    gmsh.model.setPhysicalName(2, 301, "Air")

if tag_wall:
    gmsh.model.addPhysicalGroup(2, tag_wall, 302)
    gmsh.model.setPhysicalName(2, 302, "Air_Chamber_Wall")

if tag_coolant_outside:
    gmsh.model.addPhysicalGroup(2, tag_coolant_outside, 303)
    gmsh.model.setPhysicalName(2, 303, "Coolant_Outside")

if tag_can:
    gmsh.model.addPhysicalGroup(2, tag_can, 304)
    gmsh.model.setPhysicalName(2, 304, "Al_Cassette")

for i in range(len(Fe_centres)):
    pin = i + 1
    if fuel_tags[i]:
        gmsh.model.addPhysicalGroup(2, fuel_tags[i], 100 + pin)
        gmsh.model.setPhysicalName(2, 100 + pin, f"Fuel_Pin_{pin}")
    if clad_tags[i]:
        gmsh.model.addPhysicalGroup(2, clad_tags[i], 200 + pin)
        gmsh.model.setPhysicalName(2, 200 + pin, f"Clad_Pin_{pin}")

# -------------------------------------------------
# mesh
# -------------------------------------------------
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0010)

gmsh.model.mesh.generate(2)
gmsh.write("8_rods_D5.msh")
gmsh.finalize()


