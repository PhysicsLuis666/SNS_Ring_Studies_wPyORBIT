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
from orbit.lattice import AccLattice, AccNode, AccNodeBunchTracker, AccActionsContainer
from orbit.py_linac.lattice import MarkerLinacNode, LinacTrMatricesController
import pakagesForOptimizationScripts as pk

#theta1, theta2, theta3, theta4 = sp.symbols(r"theta1, theta2, theta3, theta4") 
x,xp = sp.symbols(f"x, x'")
kx1,kx2,kx3,kx4 = sp.symbols(r"kx1, kx2, kx3, kx4")
f1,f2,f3,f4 = sp.symbols(r"f1,f2,f3,f4")
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10 = sp.symbols(r"l1, l2, l3, l4, l5, l6, l7, l8, l9, l10")

def make_drift(l):
    drift_matrix = sp.Matrix([[1, l],
                        [0, 1]])

    return drift_matrix
def make_quad_focus(k,l):
    K = k*l
    quad_matrix = sp.Matrix([[1, 0],
                             [-K,1]])
    
    return quad_matrix,
# def make_quad_defocus(f):
#     quad_matrix = sp.Matrix([[1, 0],
#                              [1/f,1]])
#     K = sp.Matrix([0,0])
#     return quad_matrix, K
def make_kick_vector(k):
    M =sp.eye(2)
    V = sp.Matrix([0,k])
    return M,V
def get_node_index(lattice, node_name):
    
    nodes = lattice.getNodes()
    i = 0
    for n in nodes:
        if n.getName() == node_name:
            return i
        
        i+=1          


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
#pk.set_kicker_strengths(lattice, "horizontal", hk_strengths)
lattice.initialize()  
 
def make_matrix(node, thin_or_thick):
    
    name = node.getName()
    typ = node.getType()
    length = node.getLength()
    
    matrix = sp.Matrix()
    if thin_or_thick == "thin":
        if typ == "drift teapot":
            matrix = [make_drift(length)]
        elif typ == "quad teapot":
            kq = node.getParam("kq")
            f = 1/(kq * length)
            matrix = [make_drift(length/2), make_quad_focus(f), make_drift(length/2)]
        elif typ == "kick teapot":
            k = 0
            if name in ["IKICKH_A10", "IKICKH_A11", "IKICKH_A12", "IKICKH_A13"]:
                if name == "IKICKH_A10":
                    k = kx1
                elif name == "IKICKH_A11":
                    k = kx2
                elif name == "IKICKH_A12":
                    k = kx3
                elif name == "IKICKH_A13":
                    k = kx4
                    
                matrix = [make_drift(length/2), make_kick_vector(k), make_drift(length/2)]
            else:
                matrix = [make_drift(length)]
        else:
            print(f"{name}: {typ}: unknown element inserting drift")
            matrix = [make_drift(length)]
    
    else:
        print("Thick elements havent been implemeneted yet")
        matrix = [sp.eye(2)]
            
    return matrix
    

        
lattice_dict = {}

for n in nodes:
    lattice_dict[n.getName()] = make_matrix(n, "thin")
    
# for n in nodes:
#     for i in lattice_dict[n.getName()]:
#         print(f"matrix in node {n.getName(), n.getType()}")
#         pprint(i)

#write function that does the track through first half of lattice
node_order = nodes[0:29]
def track_through_lattice(lattice, lattice_dict, node_order, x_i, xp_i):
    
    xxp = sp.Matrix([x_i, xp_i])
    accumulated_matrix = sp.eye(2)
    kicks_generated = []
    states_generated = []   
    
    for n in node_order:
        
        for elem in lattice_dict[n.getName()]:
            #for i in states_generated:
                #pprint(i)
            print(f"node{n.getName()}")
            pprint(elem)
            pprint(accumulated_matrix)
            if elem.shape == (2,2):
                accumulated_matrix = elem*accumulated_matrix
            elif elem.shape == (2,1):
                xxp = accumulated_matrix*xxp
                kicks_generated.append((n, accumulated_matrix, elem))
                xxp = xxp + elem
                states_generated.append(xxp)
                accumulated_matrix = sp.eye(2)
            else:
                raise RuntimeError("Unknown shape")
            
    #final apply
    xxp = accumulated_matrix*xxp
    return xxp, kicks_generated
    
x_T, xp_T = sp.symbols("xT, xpT")
print("-----------first half-----------------")
foil_xxp, array = track_through_lattice(lattice, lattice_dict, node_order, 0, 0)
    
# pprint(foil_xxp)
#print(array)
#solve for kx1, kx2
#solution = sp.solve([sp.Eq(foil_xxp[0], x_T), sp.Eq(foil_xxp[1], xp_T)], [kx1,kx2], dict=True)
solution = sp.solve([sp.Eq(foil_xxp[0], .02077959), sp.Eq(foil_xxp[1], -.0002857145)], [kx1,kx2], dict=True)

# pprint(solution)
# for i in array:
#     print(i)
# for n in nodes:
#     print(f"node {n.getName()}:{get_node_index(lattice, n.getName())}: type:{n.getType()}")

nodes_reversed_end = list(reversed(nodes))
for n in nodes_reversed_end[0:29]:
    print(f"node{n.getName()}, {n.getType()}")

print("-----------final half---------------------")

final_xxp, array2 = track_through_lattice(lattice, lattice_dict, nodes_reversed_end[0:29], 0, 0)

final_solution = sp.solve([sp.Eq(final_xxp[0],.02077959),
                           sp.Eq(final_xxp[1],-.0002857145)],
                          [kx3,kx4],
                          dict=True)#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 00:52:53 2025

@author: l5g
"""

