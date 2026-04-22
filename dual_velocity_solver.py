#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 08:13:49 2026

@author: klemenambrozic
"""

# import gmsh
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import root
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os
#import openmc
#from simulation_config import SimulationConfig
from iapws import IAPWS97



# ---------------------------------- meshing and stuff ---------------------------------------------------------------

# # Initialize gmsh
# gmsh.initialize()
# gmsh.model.add('Regular_FE')

# # Lets use opencascade cad kernel, to make the boolean operaitons easier
# occ = gmsh.model.occ



# # The [0,0] coordinate is at the centre of FE
# Fe_centres = [(-0.0255,-0.0255),
#               (-0.0255,-0.0085),
#               (-0.0255, 0.0085),
#               (-0.0255, 0.0255),
#               (-0.0085,-0.0255),
#               (-0.0085,-0.0085),
#               (-0.0085, 0.0085),
#               (-0.0085, 0.0255),
#               (0.0085,-0.0255),
#               (0.0085,-0.0085),
#               (0.0085, 0.0085),
#               (0.0085, 0.0255),
#               (0.0255,-0.0255),
#               (0.0255,-0.0085),
#               (0.0255, 0.0085),
#               (0.0255, 0.0255)]


# # --- 2. GEOMETRY PARAMETERS ---
# r_matrix =  0.0035
# r_al_clad = 0.0050

# r_clad = 0.0050

# unit_cell_size = 0.072
# box_channel_out = 0.068
# box_channel_in = 0.064 


# # --- 3. CREATE ALL OVERLAPPING PRIMITIVES ---
# # We just throw all the shapes down on top of each other.
# fluid_box = occ.addRectangle(-unit_cell_size/2, -unit_cell_size/2, 0, unit_cell_size, unit_cell_size)
# clad_out = occ.addRectangle(-box_channel_out/2,-box_channel_out/2,0,box_channel_out,box_channel_out)
# clad_in = occ.addRectangle(-box_channel_in/2,-box_channel_in/2,0,box_channel_in,box_channel_in)


# clad_disks = []
# fuel_disks = []
# for cx, cy in Fe_centres:
#     clad_disks.append(occ.addDisk(cx, cy, 0, r_al_clad, r_al_clad))
#     fuel_disks.append(occ.addDisk(cx, cy, 0, r_matrix, r_matrix))


# # --- 4. THE MAGIC FRAGMENT OPERATION ---
# # This shatters all overlapping shapes into distinct, non-overlapping pieces.
# all_tools = [(2, tag) for tag in [clad_out,clad_in]+clad_disks + fuel_disks]
# occ.fragment([(2, fluid_box)], all_tools)


# # Must synchronize before we can ask gmsh for bounding boxes
# occ.synchronize()

# # --- 5. AUTOMATICALLY SORT THE SHATTERED PIECES ---
# # Get all the new 2D surfaces created by the fragmentation
# surfaces = occ.getEntities(2)

# tag_in_water = None
# tag_out_water = None
# tag_al_box= None
# fuel_tags = {i:[] for i in range(len(Fe_centres))}
# clad_tags = {i:[] for i in range(len(Fe_centres))}



# for dim, tag in surfaces:
#     # get bounding boxes for each piece
#     xmin, ymin,zmin,xmax,ymax,zmax = gmsh.model.getBoundingBox(dim, tag)
#     width = round(xmax-xmin,5)
    
#     if width>=0.071:
#         tag_out_water = tag
#     elif width>=0.067:
#         tag_al_box = tag
#     elif width>=0.063:
#         tag_in_water = tag
#     else:
#         com_x = 0.5*(xmin+xmax)
#         com_y = 0.5*(ymin+ymax)
        
#         for i, (cx,cy) in enumerate(Fe_centres):
#             dist = math.sqrt((com_x-cx)**2 + (com_y-cy)**2)
#             if dist<=0.001:
#                 if width>=2*r_clad:
#                     clad_tags[i].append(tag)
#                 else:
#                     fuel_tags[i].append(tag)
    
    

# gmsh.model.addPhysicalGroup(2,[tag_in_water], 300, name='Inner_Coolant')
# gmsh.model.addPhysicalGroup(2,[tag_out_water], 301, name='Outer_Coolant')
# gmsh.model.addPhysicalGroup(2,[tag_al_box],400,name='Al_Cassete_Box')


# for i in range(len(Fe_centres)):
#     pin_num = i + 1
#     gmsh.model.addPhysicalGroup(2, fuel_tags[i], 100 + pin_num, name=f"Fuel_Pin_{pin_num}")
#     gmsh.model.addPhysicalGroup(2, clad_tags[i], 200 + pin_num, name=f"Clad_Pin_{pin_num}")

# # --- 7. GENERATE MESH ---
# gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0005)
# gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0001)

# print("Generating 2D Mesh...")
# gmsh.model.mesh.generate(2)
# gmsh.write("bme_assembly_4x4_fixed.msh")
# print("Saved as 'bme_assembly_4x4_fixed.msh'")

# # Visualize!
# # gmsh.fltk.run()
# # gmsh.finalize()


 


################################# fuel ltemperature and fluid flow #############################################################

class Dual_Velocity_Assembly_Solver:
    def __init__(self, config, mesh_filename):
        print("1. Loading mesh and parsing 16-pin geometry...")
        self.mesh = meshio.read(mesh_filename)
        self.nodes = self.mesh.points[:, :2]
        self.num_nodes = len(self.nodes)
        self.config = config

        self.solid_elements = []
        self.inner_fluid_elements = []
        self.outer_fluid_elements = []
        self.k_dict = {}       
        self.q_dict = {}       
        
        self.elem_pin_id = {}  # element index -> pin_id (1..16)


        cells = self.mesh.cells_dict["triangle"]
        cell_data = self.mesh.cell_data_dict["gmsh:physical"]["triangle"]
        
        for i, tag in enumerate(cell_data):
            if 101 <= tag <= 116:    # FUEL PINS
                self.solid_elements.append(i)
                self.k_dict[i] = 150.0 # If you set this to 3 (for UO2), you get a much nicer temperature profile, with Mg (150), you can barely see it 
                self.q_dict[i] = True
                self.elem_pin_id[i] = tag - 100
            elif 201 <= tag <= 216:  # CLADDING
                self.solid_elements.append(i)
                self.k_dict[i] = 200.0  
                self.q_dict[i] = False
            elif tag == 400:         # AL CASSETTE BOX
                self.solid_elements.append(i)
                self.k_dict[i] = 200.0  
                self.q_dict[i] = False
            elif tag == 300:         # INNER COOLANT
                self.inner_fluid_elements.append(i)
            elif tag == 301:         # OUTER BYPASS POOL
                self.outer_fluid_elements.append(i)

                

        # --- EXTRACT EXACT FLOW AREAS AND PERIMETERS ---
        self.A_flow_inner = 0.0
        self.A_flow_outer = 0.0
        
        for i in self.inner_fluid_elements:  # index of a triangle in the mesh
            pts = self.nodes[cells[i]]      # cells[i] gives the array of nodes IDs that triangle i is composed of, nodes[cells[i]] gives the coordinates [x, y] for each triangle ID
            self.A_flow_inner += 0.5 * abs((pts[1,0]-pts[0,0])*(pts[2,1]-pts[0,1]) - (pts[2,0]-pts[0,0])*(pts[1,1]-pts[0,1]))
            # area of the triangle using determinant

        for i in self.outer_fluid_elements:
            pts = self.nodes[cells[i]]
            self.A_flow_outer += 0.5 * abs((pts[1,0]-pts[0,0])*(pts[2,1]-pts[0,1]) - (pts[2,0]-pts[0,0])*(pts[1,1]-pts[0,1]))

        # Find all wetted boundaries by looking for edges touched by exactly 1 fluid triangle
        inner_edges = {}; outer_edges = {}
        for i in self.inner_fluid_elements:
            for j in range(3):
                edge = tuple(sorted((cells[i][j], cells[i][(j+1)%3])))
                inner_edges[edge] = inner_edges.get(edge, 0) + 1
        for i in self.outer_fluid_elements:
            for j in range(3):
                edge = tuple(sorted((cells[i][j], cells[i][(j+1)%3])))
                outer_edges[edge] = outer_edges.get(edge, 0) + 1
                
        self.P_wet_inner = sum([np.linalg.norm(self.nodes[e[0]] - self.nodes[e[1]]) for e, c in inner_edges.items() if c == 1])
        self.P_wet_outer = sum([np.linalg.norm(self.nodes[e[0]] - self.nodes[e[1]]) for e, c in outer_edges.items() if c == 1])
        
        self.D_h_inner = 4 * self.A_flow_inner / self.P_wet_inner
        self.D_h_outer = 4 * self.A_flow_outer / self.P_wet_outer

        # Identify Solid-Fluid Convective Boundaries
        self.inner_bnd_edges = []
        self.outer_bnd_edges = []
        edge_to_elements = {}
        for i, elem in enumerate(cells):   
            for j in range(3):
                edge = tuple(sorted((elem[j], elem[(j+1)%3])))
                if edge not in edge_to_elements: edge_to_elements[edge] = []
                edge_to_elements[edge].append(i)   # i will give which triangle ID is touching the edge
                
        solid_set = set(self.solid_elements)
        inner_set = set(self.inner_fluid_elements)
        outer_set = set(self.outer_fluid_elements)
        
        for edge, elems in edge_to_elements.items():
            if len(elems) == 2:
                e1, e2 = elems
                if (e1 in solid_set and e2 in inner_set) or (e2 in solid_set and e1 in inner_set):
                    self.inner_bnd_edges.append(edge)
                elif (e1 in solid_set and e2 in outer_set) or (e2 in solid_set and e1 in outer_set):
                    self.outer_bnd_edges.append(edge)

        # Track which nodes belong to which region for clean plotting
        self.inner_nodes_set = set([n for i in self.inner_fluid_elements for n in cells[i]])
        self.outer_nodes_set = set([n for i in self.outer_fluid_elements for n in cells[i]])
        self.solid_nodes_set = set([n for i in self.solid_elements for n in cells[i]])

        # --- OPTIMIZATION: PRE-ASSEMBLE THE SOLID CONDUCTION MATRIX ---
        print("2. Pre-assembling static Solid FEM Matrix (CPU Optimization)...")
        self.K_solid_base = sp.lil_matrix((self.num_nodes, self.num_nodes))
        self.F_q_unit = np.zeros(self.num_nodes) 

        

        ## check if this has to be updated to be coupled with openmc!

        for i in self.solid_elements:
            elem = cells[i]
            pts = self.nodes[elem]
            Area = 0.5 * abs((pts[1,0]-pts[0,0])*(pts[2,1]-pts[0,1]) - (pts[2,0]-pts[0,0])*(pts[1,1]-pts[0,1]))
            b = np.array([pts[1,1]-pts[2,1], pts[2,1]-pts[0,1], pts[0,1]-pts[1,1]])
            c = np.array([pts[2,0]-pts[1,0], pts[0,0]-pts[2,0], pts[1,0]-pts[0,0]])
            K_local = (self.k_dict[i] / (4 * Area)) * (np.outer(b, b) + np.outer(c, c))
            
            for m in range(3):
                for n in range(3): 
                    self.K_solid_base[elem[m], elem[n]] += K_local[m, n]
                if self.q_dict[i]: 
                    self.F_q_unit[elem[m]] += Area / 3.0
                    
        self.K_solid_base_csr = self.K_solid_base.tocsr()

        # --- 2D VISCOUS FLUID VELOCITY SHAPE ---
        self.phi_shape, self.exact_f_Re_inner, self.phi_avg = self.solve_viscous_velocity_shape()

        # TH PARAMETERS & ARRAYS
        self.L = 0.5
        self.T_inlet = config.T_inlet
        self.g = 9.81
        self.pin_list = sorted(set(self.elem_pin_id.values()))
        self.n_pins = len(self.pin_list)
        self.N_z = config.n_z

        self.dz = self.L / self.N_z
        self.z_nodes = np.linspace(0, self.L, self.N_z)
        self.r_matrix = 0.0035
        self.A_fuel = np.pi * (self.r_matrix**2)

        self.plot_T_bulk_in = np.zeros(self.N_z); self.plot_T_bulk_out = np.zeros(self.N_z)
        self.plot_v_bulk_in = np.zeros(self.N_z); self.plot_v_bulk_out = np.zeros(self.N_z)
        self.plot_T_fuel = np.zeros(self.N_z)
        self.hottest_T_nodes_2d = None
        self.hottest_v_nodes_2d = None
        self.hottest_z = 0.0

        self.qprime_pin = np.zeros((self.n_pins, self.N_z))  # W/m per pin per axial slice

        self.plot_heat_in = np.zeros(self.N_z)
        self.plot_heat_out = np.zeros(self.N_z)
        self.plot_cp_in = np.zeros(self.N_z)
        self.plot_cp_out = np.zeros(self.N_z)
        self.plot_mdot_in = 0.0
        self.plot_mdot_out = 0.0


        # blockage_ratio = 0.52
        # sigma_bottom = 1 - blockage_ratio

        self.A_holes_bottom = config.A_holes_bottom_map[os.path.basename(mesh_filename)] # Replace with the real area of the holes from your drawings [m^2]
        sigma_bottom = self.A_holes_bottom/ self.A_flow_inner

        blockage_ratio = 1 - sigma_bottom  


        K_cont = 0.5 * (1.0 - sigma_bottom) 
        K_exp = (1.0 - sigma_bottom)**2  

        K_bottom_holder = (K_cont + K_exp) #* (1.0 / sigma_bottom)**2  

        self.K_bottom = K_bottom_holder   # from Borda–Carnot or literature   2.3 
        self.K_top    = K_bottom_holder
        
        # --- PIN BOOKKEEPING ---


        print(f"Detected {self.n_pins} fuel pins: {self.pin_list}")

    def rho_old(self, T): return 1000.0 - 0.2 * (T - 20.0)
    def cp_old(self, T):  return 4184.0
    def mu_old(self, T):  return 0.001 * np.exp(-0.02*(T-20))
    def k_f_old(self, T): return 0.6    


    def rho(self, T_K, P_MPa=0.1013):
        T_K = T_K + 273.15
        # print("DEBUG rho: T_K =", T_K, "P_MPa =", P_MPa)
        T = np.asarray(T_K, dtype=float)
        rho = np.zeros_like(T, dtype=float)

        for i, Ti in np.ndenumerate(T):
            # print("DEBUG IAPWS call: T =", float(Ti), "P =", P_MPa)
            state = IAPWS97(T=float(Ti), P=float(P_MPa))
            rho[i] = state.rho
        # print("DEBUG IAPWS rho: T =", float(Ti), "cp =",  rho) 
        return rho 
    
    def mu(self, T_K, P_MPa=0.1013):
        T_K = T_K + 273.15
        T = np.asarray(T_K)
        out = np.zeros_like(T, dtype=float)
        for i, Ti in np.ndenumerate(T):
            out[i] = IAPWS97(T=Ti, P=P_MPa).mu
        # print("DEBUG IAPWS mu: T =", float(Ti), "mu =",  out)
        return out   # Pa.s

    def cp(self, T_K, P_MPa=0.1013):
        T_K = T_K + 273.15
        T = np.asarray(T_K)
        out = np.zeros_like(T, dtype=float)
        for i, Ti in np.ndenumerate(T):
            out[i] = IAPWS97(T=Ti, P=P_MPa).cp * 1e3
        # print("DEBUG IAPWS cp: T =", float(Ti), "cp =",  out) 
        return out  # J/kg/K

    def k_f(self, T_K, P_MPa=0.1013):
        T_K = T_K + 273.15
        T = np.asarray(T_K)
        out = np.zeros_like(T, dtype=float)
        for i, Ti in np.ndenumerate(T):
            out[i] = IAPWS97(T=Ti, P=P_MPa).k
        # print("DEBUG IAPWS k_f: T =", float(Ti), "k_f =",  out) 
        return out  # W/m/K

    def compute_htc(self, rho, mu, cp, k, v, D_h, L):
        Re = rho * v * D_h / mu
        Pr = cp * mu / k
        # print("Re", Re, "Pr:", Pr)
        if Re < 2300:
            Nu = 1.86 * (Re * Pr * D_h / L)**(1/3)
            Nu = max(Nu, 3.66)
        else:
            Nu = 0.023 * Re**0.8 * Pr**0.4
        # print("Nu:", Nu)
        return Nu * k / D_h, Re, Pr, Nu

    def q_vol(self, pin_id, k):
        # pin_id: 1..n_pins, k: 0..N_z-1
        return self.qprime_pin[pin_id - 1, k] / self.A_fuel   # W/m^3

    def solve_viscous_velocity_shape(self):
        print("3. Solving 2D Viscous Flow Field for exact inner friction factor...")
        K = sp.lil_matrix((self.num_nodes, self.num_nodes))
        F = np.zeros(self.num_nodes)
        node_areas = np.zeros(self.num_nodes)
        cells = self.mesh.cells_dict["triangle"]
        
        for i in self.inner_fluid_elements:
            elem = cells[i]
            pts = self.nodes[elem]
            Area = 0.5 * abs((pts[1,0]-pts[0,0])*(pts[2,1]-pts[0,1]) - (pts[2,0]-pts[0,0])*(pts[1,1]-pts[0,1]))
            b = np.array([pts[1,1]-pts[2,1], pts[2,1]-pts[0,1], pts[0,1]-pts[1,1]])
            c = np.array([pts[2,0]-pts[1,0], pts[0,0]-pts[2,0], pts[1,0]-pts[0,0]])
            K_local = (1.0 / (4 * Area)) * (np.outer(b, b) + np.outer(c, c))
            
            for m in range(3):
                node_areas[elem[m]] += Area / 3.0
                F[elem[m]] += Area / 3.0
                for n in range(3): K[elem[m], elem[n]] += K_local[m, n]

        bnd_nodes = set([n for edge in self.inner_bnd_edges for n in edge])
        for i in range(self.num_nodes):
            if i in bnd_nodes or i not in self.inner_nodes_set:
                K[i, :] = 0; K[i, i] = 1.0; F[i] = 0.0
                
        phi = spla.spsolve(K.tocsr(), F)
        phi_avg = np.sum(phi * node_areas) / self.A_flow_inner
        exact_f_Re = (2 * self.D_h_inner**2) / phi_avg
        print(f"   -> Exact Inner Channel f*Re: {exact_f_Re:.2f}")
        return phi, exact_f_Re, phi_avg

    def solve_2d_slice(self, z, T_inner, h_inner, T_outer, h_outer):
        cells = self.mesh.cells_dict["triangle"]
        
        K_conv = sp.lil_matrix((self.num_nodes, self.num_nodes))
        
        #F = self.F_q_unit * self.q_vol(z)
        k = min(int(z / self.dz), self.N_z - 1)
        F = np.zeros(self.num_nodes)

        for i in self.solid_elements:
            if not self.q_dict[i]:
                continue
            pin = self.elem_pin_id[i]          # 1..n_pins
            qppp = self.q_vol(pin, k)          # W/m^3

            elem = cells[i]
            pts = self.nodes[elem]
            Area = 0.5 * abs((pts[1,0]-pts[0,0])*(pts[2,1]-pts[0,1]) - (pts[2,0]-pts[0,0])*(pts[1,1]-pts[0,1]))

            for m in range(3):
                F[elem[m]] += qppp * Area / 3.0

        

        # Inner Coolant Convection
        for n1, n2 in self.inner_bnd_edges:
            L_edge = np.linalg.norm(self.nodes[n1] - self.nodes[n2])
            conv_mat = (h_inner * L_edge / 6.0) * np.array([[2, 1], [1, 2]])
            K_conv[n1, n1] += conv_mat[0,0]; K_conv[n1, n2] += conv_mat[0,1]
            K_conv[n2, n1] += conv_mat[1,0]; K_conv[n2, n2] += conv_mat[1,1]
            F[n1] += (h_inner * T_inner * L_edge / 2.0); F[n2] += (h_inner * T_inner * L_edge / 2.0)

        # Outer Pool Convection
        for n1, n2 in self.outer_bnd_edges:
            L_edge = np.linalg.norm(self.nodes[n1] - self.nodes[n2])
            conv_mat = (h_outer * L_edge / 6.0) * np.array([[2, 1], [1, 2]])
            K_conv[n1, n1] += conv_mat[0,0]; K_conv[n1, n2] += conv_mat[0,1]
            K_conv[n2, n1] += conv_mat[1,0]; K_conv[n2, n2] += conv_mat[1,1]
            F[n1] += (h_outer * T_outer * L_edge / 2.0); F[n2] += (h_outer * T_outer * L_edge / 2.0)

        # Fix unused nodes for clean plotting
        for i in range(self.num_nodes):
            if i not in self.solid_nodes_set:
                K_conv[i, i] = 1.0
                F[i] = T_inner if i in self.inner_nodes_set else T_outer

        K_final = self.K_solid_base_csr + K_conv.tocsr()
        T_nodes = spla.spsolve(K_final, F)
        
        heat_to_inner = sum([h_inner * np.linalg.norm(self.nodes[n1] - self.nodes[n2]) * (((T_nodes[n1] + T_nodes[n2])/2.0) - T_inner) for n1, n2 in self.inner_bnd_edges])
        heat_to_outer = sum([h_outer * np.linalg.norm(self.nodes[n1] - self.nodes[n2]) * (((T_nodes[n1] + T_nodes[n2])/2.0) - T_outer) for n1, n2 in self.outer_bnd_edges])

        outer_wall_avg = np.mean([
            (T_nodes[n1] + T_nodes[n2]) / 2.0
            for n1, n2 in self.outer_bnd_edges
        ])

        inner_wall_avg = np.mean([
            (T_nodes[n1] + T_nodes[n2]) / 2.0
            for n1, n2 in self.inner_bnd_edges
        ])

        print(
            f"z={z:.5f}, "
            f"T_inner={T_inner:.6f}, T_outer={T_outer:.6f}, "
            f"Twall_inner_avg={outer_wall_avg:.6f}, Twall_outer_avg={outer_wall_avg:.6f}, "
            f"heat_in={heat_to_inner:.6e}, heat_out={heat_to_outer:.6e}",
            flush=True
        )

            
        return heat_to_inner, heat_to_outer, np.max(T_nodes), T_nodes

    # def evaluate_flow(self, v_guesses, save_results=False):
    #     v_in, v_out = v_guesses[0], v_guesses[1]
        
    #     # Heavily penalize negative or zero velocities to keep the solver on track
    #     if v_in <= 1e-6 or v_out <= 1e-6: 
    #         return [1e6, 1e6] 
        

            
    #     # Hydrostatic pressure setup
    #     rho_avg = 1000.0          # kg/m³ ≈ liquid water
    #     P_surface = 0.1013        # MPa at pool surface
    #     g = 9.81                  # m/s²

    #     # z = 0 at core top → depth from pool surface: h = 5.0 + z
    #     # z_nodes are in meters from core top downward
    #     # p_MPa(z) = P_surface + (rho_avg * g * (5.0 + z)) * 1e-6
    #     P_z = P_surface + rho_avg * g * (5.0 + self.z_nodes) * 1e-6   # array of MPa

    #     # Inlet properties at **core‑top** pressure
    #     p_inlet_MPa = P_z[0]   # pressure at core top (z=0)
    #     rho_inlet = self.rho(self.T_inlet, P_MPa=p_inlet_MPa)  # T in °C, function 
        
    #     T_bulk_in = self.T_inlet; T_bulk_out = self.T_inlet
    #     # T_bulk_in = np.clip(T_bulk_in, 273.15, 800.0)
    #     # T_bulk_out = np.clip(T_bulk_out, 273.15, 800.0)

    #     rho_inlet = self.rho(self.T_inlet)
    #     m_dot_in = rho_inlet * self.A_flow_inner * v_in
    #     m_dot_out = rho_inlet * self.A_flow_outer * v_out
        
    #     buoyancy_in, friction_in = 0.0, 0.5 * 0.5 * rho_inlet * v_in**2
    #     buoyancy_out, friction_out = 0.0, 0.5 * 0.5 * rho_inlet * v_out**2

        
    #     for i, z in enumerate(self.z_nodes):
    #         p_MPa = P_z[i]   # local pressure in MPa at current z

    #         # T_bulk_in, T_bulk_out are in °C
    #         rho_in_l = self.rho(T_bulk_in, P_MPa=p_MPa)
    #         mu_in_l  = self.mu( T_bulk_in, P_MPa=p_MPa)
    #         cp_in_l  = self.cp( T_bulk_in, P_MPa=p_MPa)

    #         rho_out_l = self.rho(T_bulk_out, P_MPa=p_MPa)
    #         mu_out_l  = self.mu( T_bulk_out, P_MPa=p_MPa)
    #         cp_out_l  = self.cp( T_bulk_out, P_MPa=p_MPa)

    #         # # Fluid Properties
    #         # rho_in_l, mu_in_l, cp_in_l = self.rho(T_bulk_in), self.mu(T_bulk_in), self.cp(T_bulk_in)
    #         # rho_out_l, mu_out_l, cp_out_l = self.rho(T_bulk_out), self.mu(T_bulk_out), self.cp(T_bulk_out)
            
    #         v_in_local = m_dot_in / (rho_in_l * self.A_flow_inner)
    #         v_out_local = m_dot_out / (rho_out_l * self.A_flow_outer)
            
    #         # Friction
    #         Re_in = rho_in_l * v_in_local * self.D_h_inner / mu_in_l
    #         f_in = self.exact_f_Re_inner / max(Re_in, 1) if Re_in < 2000 else 0.316 * Re_in**(-0.25)
    #         friction_in += f_in * (self.dz / self.D_h_inner) * 0.5 * rho_in_l * v_in_local**2
            
    #         Re_out = rho_out_l * v_out_local * self.D_h_outer / mu_out_l
    #         f_out = 96.0 / max(Re_out, 1) if Re_out < 2000 else 0.316 * Re_out**(-0.25)
    #         friction_out += f_out * (self.dz / self.D_h_outer) * 0.5 * rho_out_l * v_out_local**2
            
    #         # Convection & Solid Solve
    #         h_in = 4.36 * self.k_f(T_bulk_in, P_MPa=p_MPa) / self.D_h_inner 
    #         h_out = 4.36 * self.k_f( T_bulk_out, P_MPa=p_MPa) / self.D_h_outer 
    #         # h_in, Re_in, Pr_in, Nu_in = self.compute_htc(
    #         #     rho_in_l, mu_in_l, cp_in_l, self.k_f(T_bulk_in),
    #         #     v_in_local, self.D_h_inner, self.L
    #         # )

    #         # h_out, Re_out, Pr_out, Nu_out = self.compute_htc(
    #         #     rho_out_l, mu_out_l, cp_out_l, self.k_f(T_bulk_out),
    #         #     v_out_local, self.D_h_outer, self.L
    #         # )
                        
    #         heat_in, heat_out, T_max, T_nodes_2d = self.solve_2d_slice(z, T_bulk_in, h_in, T_bulk_out, h_out)
            
    #         if save_results:
    #             self.plot_T_bulk_in[i] = T_bulk_in; self.plot_T_bulk_out[i] = T_bulk_out
    #             self.plot_v_bulk_in[i] = v_in_local; self.plot_v_bulk_out[i] = v_out_local
    #             self.plot_T_fuel[i] = T_max
    #             if i == self.N_z // 2:
    #                 self.hottest_T_nodes_2d = T_nodes_2d
    #                 self.hottest_z = z
    #                 self.hottest_v_nodes_2d = v_in_local * (self.phi_shape / self.phi_avg)
            
    #         # Buoyancy & Energy Equations
    #         buoyancy_in += self.g * (rho_inlet - rho_in_l) * self.dz
    #         buoyancy_out += self.g * (rho_inlet - rho_out_l) * self.dz
            
    #         T_bulk_in += (heat_in * self.dz) / (m_dot_in * cp_in_l)
    #         T_bulk_out += (heat_out * self.dz) / (m_dot_out * cp_out_l)

    #         dT_in = (heat_in * self.dz) / (m_dot_in * cp_in_l)
    #         dT_out = (heat_out * self.dz) / (m_dot_out * cp_out_l)

    #         print(
    #             f"z={z:.5f}, "
    #             f"heat_in={heat_in:.6e}, heat_out={heat_out:.6e}, "
    #             f"m_dot_in={m_dot_in:.6e}, m_dot_out={m_dot_out:.6e}, "
    #             f"dT_in={dT_in:.6e}, dT_out={dT_out:.6e}",
    #             flush=True
    #         )
    #     friction_in += 1.0 * 0.5 * rho_in_l * v_in_local**2
    #     friction_out += 1.0 * 0.5 * rho_out_l * v_out_local**2
        
    #     # We must return an array of the two residuals.
    #     print(f'Convergence - Inside flow: {buoyancy_in - friction_in}, Outside flow: {buoyancy_out - friction_out}')
    #     return [buoyancy_in - friction_in, buoyancy_out - friction_out]



    def evaluate_flow(self, v_guesses, save_results=False):
        v_in, v_out = v_guesses[0], v_guesses[1]

        if v_in <= 1e-6 or v_out <= 1e-6:
            return [1e6, 1e6]

        # Hydrostatic pressure profile
        rho_avg = 1000.0      # kg/m^3
        P_surface = 0.1013    # MPa
        P_z = P_surface + rho_avg * self.g * (5.0 + (self.L - self.z_nodes)) * 1e-6

        p_inlet_MPa = P_z[0]

        T_bulk_in = self.T_inlet
        T_bulk_out = self.T_inlet

        rho_inlet = float(self.rho(self.T_inlet, P_MPa=p_inlet_MPa))
        m_dot_in = rho_inlet * self.A_flow_inner * v_in
        m_dot_out = rho_inlet * self.A_flow_outer * v_out


        buoyancy_in = 0.0

        # inlet kinetic + bottom grid loss
        friction_in = 0.5 * rho_inlet * v_in**2 * (0.5 + self.K_bottom)

        # buoyancy_in, friction_in = 0.0, 0.25 * rho_inlet * v_in**2
        buoyancy_out, friction_out = 0.0, 0.25 * rho_inlet * v_out**2

        for i, z in enumerate(self.z_nodes):
            p_MPa = P_z[i]

            rho_in_l = float(self.rho(T_bulk_in, P_MPa=p_MPa))
            mu_in_l  = float(self.mu(T_bulk_in, P_MPa=p_MPa))
            cp_in_l  = float(self.cp(T_bulk_in, P_MPa=p_MPa))

            rho_out_l = float(self.rho(T_bulk_out, P_MPa=p_MPa))
            mu_out_l  = float(self.mu(T_bulk_out, P_MPa=p_MPa))
            cp_out_l  = float(self.cp(T_bulk_out, P_MPa=p_MPa))

            v_in_local = m_dot_in / (rho_in_l * self.A_flow_inner)
            v_out_local = m_dot_out / (rho_out_l * self.A_flow_outer)

            Re_in = rho_in_l * v_in_local * self.D_h_inner / mu_in_l

            f_in = self.exact_f_Re_inner / max(Re_in, 1.0) if Re_in < 2000 else 0.316 * Re_in**(-0.25)
            friction_in += f_in * (self.dz / self.D_h_inner) * 0.5 * rho_in_l * v_in_local**2
        

            Re_out = rho_out_l * v_out_local * self.D_h_outer / mu_out_l
            f_out = 96.0 / max(Re_out, 1.0) if Re_out < 2000 else 0.316 * Re_out**(-0.25)
            friction_out += f_out * (self.dz / self.D_h_outer) * 0.5 * rho_out_l * v_out_local**2

            h_in = 4.36 * float(self.k_f(T_bulk_in, P_MPa=p_MPa)) / self.D_h_inner
            h_out = 4.36 * float(self.k_f(T_bulk_out, P_MPa=p_MPa)) / self.D_h_outer

            heat_in, heat_out, T_max, T_nodes_2d = self.solve_2d_slice(
                z, T_bulk_in, h_in, T_bulk_out, h_out
            )

            if save_results:
                self.plot_T_bulk_in[i] = T_bulk_in
                self.plot_T_bulk_out[i] = T_bulk_out
                self.plot_v_bulk_in[i] = v_in_local
                self.plot_v_bulk_out[i] = v_out_local
                self.plot_T_fuel[i] = T_max

                self.plot_heat_in[i] = heat_in
                self.plot_heat_out[i] = heat_out
                self.plot_cp_in[i] = cp_in_l
                self.plot_cp_out[i] = cp_out_l
                self.plot_mdot_in = m_dot_in
                self.plot_mdot_out = m_dot_out


                if i == self.N_z // 2:
                    self.hottest_T_nodes_2d = T_nodes_2d
                    self.hottest_z = z
                    self.hottest_v_nodes_2d = v_in_local * (self.phi_shape / self.phi_avg)
                    

            buoyancy_in += self.g * (rho_inlet - rho_in_l) * self.dz
            buoyancy_out += self.g * (rho_inlet - rho_out_l) * self.dz

            dT_in = (heat_in * self.dz) / (m_dot_in * cp_in_l)
            dT_out = (heat_out * self.dz) / (m_dot_out * cp_out_l)

            print(
                f"z={z:.5f}, "
                f"heat_in={heat_in:.6e}, heat_out={heat_out:.6e}, "
                f"m_dot_in={m_dot_in:.6e}, m_dot_out={m_dot_out:.6e}, "
                f"dT_in={dT_in:.6e}, dT_out={dT_out:.6e}",
                flush=True
            )
            print("sum qprime over pins * dz =", np.sum(self.qprime_pin[:, i]) * self.dz)

            T_bulk_in += dT_in

            # freeze outer coolant for stability
            T_bulk_out += dT_out

        friction_in += 0.5 * rho_in_l * v_in_local**2 * (1.0 + self.K_top)
        friction_out += 0.5 * rho_out_l * v_out_local**2


        print(f"Convergence - Inside flow: {buoyancy_in - friction_in}, Outside flow: {buoyancy_out - friction_out}")
        return [buoyancy_in - friction_in, buoyancy_out - friction_out]

    def set_power_profile(self, qprime_pin):
        """
        qprime_pin: array-like of shape (n_pins, N_z), units W/m
        """
        qprime_pin = np.asarray(qprime_pin, dtype=float)
        if qprime_pin.shape != (self.n_pins, self.N_z):
            raise ValueError(f"Expected shape {(self.n_pins, self.N_z)} but got {qprime_pin.shape}")
        self.qprime_pin[:, :] = qprime_pin


    def assign_qprime_by_centers(self, qprime_values, rod_centers):
        """
        Match OpenMC qprime values to mesh pins by closest rod center distance.
        
        Args:
            qprime_values: (n_pins, N_z) array [W/m] from OpenMC tallies
            rod_centers: (n_pins, 2) array of (x,y) rod centers from rod_mapping
        """
        if qprime_values.shape[0] != len(rod_centers):
            raise ValueError(f"Pin count mismatch: {qprime_values.shape[0]} vs {len(rod_centers)}")
        
        # Compute pin center positions from your FEM mesh
        pin_mesh_centers = self.compute_pin_centers()  # Implement below

        rod_centers = rod_centers * 1e-2
        xmin = rod_centers[:, 0].min()
        xmax = rod_centers[:, 0].max()

        ymin = rod_centers[:, 1].min()
        ymax = rod_centers[:, 1].max()

        center_x = 0.5 * (xmin + xmax)
        center_y = 0.5 * (ymin + ymax)

        assy_center = np.array([center_x, center_y])
        rod_centers = rod_centers - assy_center
                
        # Match each mesh pin to closest OpenMC rod by Euclidean distance
        distances = np.linalg.norm(pin_mesh_centers[:, None] - rod_centers[None, :], axis=2)
        pin_to_rod_idx = np.argmin(distances, axis=1)
        
        # Assign qprime values by nearest neighbor matching
        matched_qprime = qprime_values[pin_to_rod_idx, :]
        
        self.set_power_profile(matched_qprime)
        print(f"Assigned qprime: {matched_qprime.shape}, max={np.max(matched_qprime):.1f} W/m")

        
        print("pin_to_rod_idx =", pin_to_rod_idx)
        print("unique matches =", np.unique(pin_to_rod_idx))
        print("n unique =", len(np.unique(pin_to_rod_idx)), "out of", len(pin_to_rod_idx))


    def compute_pin_centers(self):
        pin_centers = {}

        for i in self.solid_elements:
            if i not in self.elem_pin_id:
                continue

            pin = self.elem_pin_id[i]
            elem = self.mesh.cells_dict["triangle"][i]
            pts = self.nodes[elem]

            center = np.mean(pts, axis=0)
            area = 0.5 * abs(
                (pts[1,0]-pts[0,0])*(pts[2,1]-pts[0,1]) -
                (pts[2,0]-pts[0,0])*(pts[1,1]-pts[0,1])
            )

            if pin not in pin_centers:
                pin_centers[pin] = {'sum': np.zeros(2), 'area': 0.0}

            pin_centers[pin]['sum'] += center * area
            pin_centers[pin]['area'] += area

        centers_array = np.empty((self.n_pins, 2))
        for j, pin in enumerate(self.pin_list):
            centers_array[j] = pin_centers[pin]['sum'] / pin_centers[pin]['area']


        return centers_array

    def plot_results(self, out_dir=None, case_name="case"):
        # default: save where the statepoint/qprime live
        if out_dir is None:
            out_dir = getattr(self, "out_dir", os.getcwd())

        os.makedirs(out_dir, exist_ok=True)
        fname = os.path.join(out_dir, f"Solver_Results_{case_name}.png")

        fig = plt.figure(figsize=(16, 12))
        triang = mtri.Triangulation(
            self.nodes[:, 0] * 1000,
            self.nodes[:, 1] * 1000,
            self.mesh.cells_dict["triangle"]
        )

        ax1 = fig.add_subplot(221)
        ax1.plot(self.z_nodes, self.plot_T_bulk_in, label='Inner Coolant Temp', linewidth=2)
        ax1.plot(self.z_nodes, self.plot_T_bulk_out, label='Outer Bypass Temp', linewidth=2, linestyle='--')
        ax1.plot(self.z_nodes, self.plot_T_fuel, label='Max Fuel Temp', linewidth=2)
        ax1.set_title('Coupled Axial Temperatures'); ax1.grid(True); ax1.legend()

        ax2 = fig.add_subplot(222)
        ax2.plot(self.z_nodes, self.plot_v_bulk_in, label='Inner Cassette Velocity', linewidth=2)
        ax2.plot(self.z_nodes, self.plot_v_bulk_out, label='Outer Bypass Velocity', linewidth=2, linestyle='--')
        ax2.set_title('Dual Fluid Acceleration'); ax2.grid(True); ax2.legend()

        ax3 = fig.add_subplot(223)
        tpc_t = ax3.tripcolor(triang, self.hottest_T_nodes_2d, shading='flat', cmap='inferno')
        ax3.set_aspect('equal'); ax3.set_title(f'Solid Temp Field at z={self.hottest_z:.2f}m')
        fig.colorbar(tpc_t, ax=ax3, label='Temperature [°C]')

        ax4 = fig.add_subplot(224)
        tpc_v = ax4.tripcolor(triang, self.hottest_v_nodes_2d, shading='flat', cmap='viridis')
        ax4.set_aspect('equal'); ax4.set_title(f'Inner Viscous Profile at z={self.hottest_z:.2f}m')
        fig.colorbar(tpc_v, ax=ax4, label='Velocity [m/s]')

        plt.tight_layout()
        fig.savefig(fname, dpi=300)     # <-- save deterministically
        plt.close(fig)                  # <-- important on cluster
        print(f"Saved figure: {fname}")

    def run_old(self, case_name="case", out_dir=None, initial_guess=None):
        if out_dir is None:
            out_dir = getattr(self, "out_dir", os.getcwd())
        self.out_dir = out_dir

        if initial_guess is None:
            P = getattr(self, "total_power", None)
            if P is None:
                dz = self.L / self.N_z
                P = float(np.sum(self.qprime_pin) * dz) if hasattr(self, "qprime_pin") else self.config.P_core_W

            P_ref = getattr(self.config, "P_ref_W", 100_000.0)
            scale = np.sqrt(max(P, 1.0) / P_ref)

            v_in_ref  = getattr(self.config, "v_in_ref", 0.1)
            v_out_ref = getattr(self.config, "v_out_ref", 0.05)
            v_min = getattr(self.config, "v_min", 1e-6)

            initial_guess = [
                max(v_in_ref * scale, 10 * v_min),
                max(v_out_ref * scale, 10 * v_min)
            ]

        guess_list = [
            initial_guess,
            [0.1, 0.05],
            [0.05, 0.02],
            [0.02, 0.01],
        ]

        sol = None
        last_trial = None
        for guess in guess_list:
            print(f"Trying initial guess {guess} for {case_name}")
            trial = root(self.evaluate_flow, guess, method="hybr")
            last_trial = trial
            if trial.success:
                sol = trial
                break

        if sol is None:
            sol = last_trial

        if sol.success:
            print(f"\n✅ CONVERGED! {case_name}")
            print(f"Inner Cassette Velocity: {sol.x[0]:.6f} m/s")
            print(f"Outer Bypass Velocity:   {sol.x[1]:.6f} m/s")
            self.evaluate_flow(sol.x, save_results=True)
            self.plot_results(out_dir=out_dir, case_name=case_name)
        else:
            print(f"\n❌ Solver failed to converge ({case_name}).")
            print("Initial guess was:", initial_guess)
            print(sol.message)

        return sol



    def run(self, case_name="case", out_dir=None, initial_guess=None):
        if out_dir is None:
            out_dir = getattr(self, "out_dir", os.getcwd())
        self.out_dir = out_dir

        if initial_guess is None:
            P = getattr(self, "total_power", None)
            if P is None:
                dz = self.L / self.N_z
                P = float(np.sum(self.qprime_pin) * dz) if hasattr(self, "qprime_pin") else self.config.P_core_W

            P_ref = getattr(self.config, "P_ref_W", 100_000.0)
            scale = np.sqrt(max(P, 1.0) / P_ref)

            v_in_ref  = getattr(self.config, "v_in_ref", 0.1)
            v_out_ref = getattr(self.config, "v_out_ref", 0.05)
            v_min = getattr(self.config, "v_min", 1e-6)

            initial_guess = [
                max(v_in_ref * scale, 10 * v_min),
                max(v_out_ref * scale, 10 * v_min)
            ]
            
        guess_list = [
            initial_guess,                # physics-based
            [0.1, 0.05],
            [0.02, 0.005],                # moderate fallback
            [0.01, 0.003],                # safe lower bound
        ]
        tol_res = getattr(self.config, "tol_flow_residual", 1e-6)

        best_trial = None
        best_residual = np.inf

        for guess in guess_list:
            print(f"Trying initial guess {guess} for {case_name}")

            trial = root(self.evaluate_flow, guess, method="hybr")

            # evaluate residual at returned point
            res = np.asarray(self.evaluate_flow(trial.x, save_results=False), dtype=float)
            res_norm = np.linalg.norm(res, ord=np.inf)

            print(f"   success={trial.success}, residual_inf={res_norm:.6e}")
            print(f"   x = {trial.x}")
            print(f"   message = {trial.message}")

            if res_norm < best_residual:
                best_residual = res_norm
                best_trial = trial

            # accept either official scipy success OR small enough residual
            if trial.success or res_norm < tol_res:
                best_trial = trial
                best_residual = res_norm
                break

        sol = best_trial

        if sol is None:
            raise RuntimeError(f"No solver trial produced a result for {case_name}")

        # final acceptance based on residual, not only scipy success
        final_res = np.asarray(self.evaluate_flow(sol.x, save_results=False), dtype=float)
        final_res_norm = np.linalg.norm(final_res, ord=np.inf)

        if final_res_norm < tol_res:
            print(f"\n✅ CONVERGED (residual-based)! {case_name}")
        elif sol.success:
            print(f"\n⚠️ CONVERGED (scipy flag only)! {case_name}")
        else:
            print(f"\n❌ Solver failed ({case_name})")

        print(f"Inner Cassette Velocity: {sol.x[0]:.6f} m/s")
        print(f"Outer Bypass Velocity:   {sol.x[1]:.6f} m/s")
        print(f"Final residual inf-norm: {final_res_norm:.6e}")
        print(f"Solver message: {sol.message}")

        if final_res_norm < tol_res or sol.success:
            self.evaluate_flow(sol.x, save_results=True)
            self.plot_results(out_dir=out_dir, case_name=case_name)

        return sol