#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 21:29:30 2025

@author: l5g
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 15:11:07 2025

@author: l5g
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 15:19:58 2025

@author: l5g
"""
#this script for now just gets the matricies straight from PyORBIT
import os 
import sys
import time as t
import numpy as np
import matplotlib.pyplot as plt
import scipy as spy

from sympy import init_printing, pprint
init_printing(use_unicode=True, wrap_line=False)
from sympy import symbols
from sympy.matrices import Matrix
from sympy import eye, simplify, MatrixSymbol
import sympy as sp
from numpy import sign

sys.path.insert(0, "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/modules_to_use")
from orbit.core.bunch import Bunch, BunchTwissAnalysis
from orbit.teapot.teapot_matrix_lattice import TEAPOT_MATRIX_Lattice, TEAPOT_Lattice
from orbit.matrix_lattice import BaseMATRIX, MATRIX_Lattice
from orbit.teapot_base import MatrixGenerator
from orbit.lattice import AccLattice, AccNode, AccNodeBunchTracker, AccActionsContainer
from orbit.py_linac.lattice import MarkerLinacNode, LinacTrMatricesController
import pakagesForOptimizationScripts as pk

#theta1, theta2, theta3, theta4 = sp.symbols(r"theta1, theta2, theta3, theta4") 
x,xp = sp.symbols(f"x, x'")
kx1,kx2,kx3,kx4 = sp.symbols(r"kh1, kh2, kh3, kh4")
ky1,ky2,ky3,ky4 = sp.symbols(r'kv1, kv2, kv3, kv4')



def make_kick_vector(k):
    
    return sp.Matrix([0,k])
def get_node_index(lattice, node_name):
    
    nodes = lattice.getNodes()
    i = 0
    for n in nodes:
        if n.getName() == node_name:
            return i
        
        i+=1          


def get_matrix(matrix):
    
    n = matrix.size()[0]
    m = matrix.size()[1]
    
    matrix_test = sp.zeros(n,m)
    for i in range(n):
        for j in range(m):
            matrix_test[i,j] = matrix.get(i,j)
    return matrix_test
#pprint(kicker_vector)
bunch1 = Bunch()
pk.do_bunch(bunch1, 1.3)
lattice = pk.build_injection_region_lattice()
nodes = lattice.getNodes()
# for n in nodes:
#     print(f"{n.getName()}: {n.getType()}")
hkicker_nodes = pk.get_kickers(lattice, "horizontal")
hkicks = pk.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
hk_strengths = pk.obtain_scaled_kicker_strengths(lattice, "horizontal", hkicks, [0,0])
pk.set_kicker_strengths(lattice, "horizontal", hk_strengths)
lattice.initialize()  
matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch1)
matrix_lattice.initialize()
matrix_generator = MatrixGenerator()
matrix_nodes = matrix_lattice.getNodes()




for matrixNodes in matrix_nodes:
    print(f"node :{matrixNodes.getName()}, index {get_node_index(matrix_lattice, matrixNodes.getName())}")
    matrix = matrixNodes.getMatrix()
    trf_matrix = get_matrix(matrix)
    pprint(trf_matrix)

print(dir(matrix_nodes))
kickh_a10 = matrix_lattice.getNodeForName("IKICKH_A10")

lattice_dict = {}
print(hk_strengths)
#for mtrxNodes in matrix_nodes:
    
    

import sympy as sp

kx1, kx2, kx3, kx4 = sp.symbols("kx1 kx2 kx3 kx4")

# Map kicker names to symbols
kicker_symbols = {
    "IKICKH_A10": kx1,
    "IKICKH_A11": kx2,
    "IKICKH_A12": kx3,
    "IKICKH_A13": kx4
}

lattice_dict = {}

for node in matrix_nodes:
    name = node.getName()
    matrix = node.getMatrix()
    n, m = matrix.size()
    M_sym = sp.zeros(n, m)
    
    for i in range(n):
        for j in range(m):
            val = matrix.get(i, j)
            
            # 🔹 Only replace kick values if node is a horizontal kicker
            if name in kicker_symbols:
                # Detect the kick entry — typically (0,1) or (1,0) depending on how PyORBIT encodes it
                # Most often, the kick term appears in the (1,0) position in thin-lens form
                if abs(val) > 0 and i == 1 and j == 6:
                    M_sym[i, j] = kicker_symbols[name]
                else:
                    M_sym[i, j] = val
            else:
                M_sym[i, j] = val
    
    lattice_dict[name] = [M_sym]






# def make_matrix(node, thin_or_thick):
    
#     name = node.getName()
#     typ = node.getType()
#     length = node.getLength()
    
#     matrix = sp.Matrix()
#     if thin_or_thick == "thin":
#         if typ == "drift teapot":
#             matrix = [make_drift(length)]
#         elif typ == "quad teapot":
#             kq = node.getParam("kq")
#             f = 1/(kq * length)
#             matrix = [make_drift(length/2), make_quad_focus(f), make_drift(length/2)]
#         elif typ == "kick teapot":
#             k = 0
#             if name in ["IKICKH_A10", "IKICKH_A11", "IKICKH_A12", "IKICKH_A13"]:
#                 if name == "IKICKH_A10":
#                     k = kx1
#                 elif name == "IKICKH_A11":
#                     k = kx2
#                 elif name == "IKICKH_A12":
#                     k = kx3
#                 elif name == "IKICKH_A13":
#                     k = kx4
                    
#                 matrix = [make_drift(length/2), make_kick_vector(k), make_drift(length/2)]
#             else:
#                 matrix = [make_drift(length)]
#         else:
#             print(f"{name}: {typ}: unknown element inserting drift")
#             matrix = [make_drift(length)]
    
#     else:
#         print("Thick elements havent been implemeneted yet")
#         matrix = [sp.eye(2)]
            
#     return matrix
    

        
# lattice_dict = {}

# for n in nodes:
#     lattice_dict[n.getName()] = make_matrix(n, "thin")
    
# # for n in nodes:
# #     for i in lattice_dict[n.getName()]:
# #         print(f"matrix in node {n.getName(), n.getType()}")
# #         pprint(i)

# #write function that does the track through first half of lattice
# node_order = nodes[0:29]
# def track_through_lattice(lattice, lattice_dict, node_order, x_i, xp_i):
    
#     xxp = sp.Matrix([x_i, xp_i])
#     accumulated_matrix = sp.eye(2)
#     kicks_generated = []
#     states_generated = []   
    
#     for n in node_order:
        
#         for elem in lattice_dict[n.getName()]:
#             #for i in states_generated:
#                 #pprint(i)
#             print(f"node{n.getName()}")
#             pprint(elem)
#             pprint(accumulated_matrix)
#             if elem.shape == (2,2):
#                 accumulated_matrix = elem*accumulated_matrix
#             elif elem.shape == (2,1):
#                 xxp = accumulated_matrix*xxp
#                 kicks_generated.append((n, accumulated_matrix, elem))
#                 xxp = xxp + elem
#                 states_generated.append(xxp)
#                 accumulated_matrix = sp.eye(2)
#             else:
#                 raise RuntimeError("Unknown shape")
            
#     #final apply
#     xxp = accumulated_matrix*xxp
#     return xxp, kicks_generated
    
# x_T, xp_T = sp.symbols("xT, xpT")
# print("-----------first half-----------------")
# foil_xxp, array = track_through_lattice(lattice, lattice_dict, node_order, 0, 0)
    
# # pprint(foil_xxp)
# #print(array)
# #solve for kx1, kx2
# #solution = sp.solve([sp.Eq(foil_xxp[0], x_T), sp.Eq(foil_xxp[1], xp_T)], [kx1,kx2], dict=True)
# solution = sp.solve([sp.Eq(foil_xxp[0], .02077959), sp.Eq(foil_xxp[1], -.0002857145)], [kx1,kx2], dict=True)

# # pprint(solution)
# # for i in array:
# #     print(i)
# # for n in nodes:
# #     print(f"node {n.getName()}:{get_node_index(lattice, n.getName())}: type:{n.getType()}")

# nodes_reversed_end = list(reversed(nodes))
# for n in nodes_reversed_end[0:29]:
#     print(f"node{n.getName()}, {n.getType()}")

# print("-----------final half---------------------")

# final_xxp, array2 = track_through_lattice(lattice, lattice_dict, nodes_reversed_end[0:29], 0, 0)

# final_solution = sp.solve([sp.Eq(final_xxp[0],.02077959),
#                            sp.Eq(final_xxp[1],-.0002857145)],
#                           [kx3,kx4],
#                           dict=True)