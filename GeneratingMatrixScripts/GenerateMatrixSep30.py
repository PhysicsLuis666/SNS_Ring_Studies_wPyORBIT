#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 14:10:04 2025

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
ky1,ky2,ky3,ky4 = sp.symbols(r"ky1,ky2,ky3,ky4")
f1,f2,f3,f4 = sp.symbols(r"f1,f2,f3,f4")
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10 = sp.symbols(r"l1, l2, l3, l4, l5, l6, l7, l8, l9, l10")
sp.pprint(kx1)

bunch1 = Bunch()
pk.do_bunch(bunch1, 1.3)
lattice = pk.build_injection_region_lattice()
nodes = lattice.getNodes()
hkicker_nodes = pk.get_kickers(lattice, "horizontal")
hkicks = pk.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
hk_strengths = pk.obtain_scaled_kicker_strengths(lattice, "horizontal", hkicks, [0,0])
pk.set_kicker_strengths(lattice, "horizontal", hk_strengths)


def make_drift(l):
    drift_matrix = sp.Matrix([[1, l],
                        [0, 1]])
    return drift_matrix
def make_quad_focus(f):
    quad_matrix = sp.Matrix([[1, 0],
                             [-1/f,1]])
    return quad_matrix
def make_quad_defocus(f):
    quad_matrix = sp.Matrix([[1, 0],
                             [1/f,1]])
    return quad_matrix
def make_kick_vector(k):
    
    if k != 0:
        vector = sp.Matrix([0, k])
    else:
        vector = sp.eye(2)
    return vector
def get_node_index(lattice, node_name):
    
    nodes = lattice.getNodes()
    i = 0
    for n in nodes:
        if n.getName() == node_name:
            return i
        i+=1

for n in nodes:
    print(f"{n.getName()}:{n.getType()}, length:{n.getLength()}")
    
def make_matrix(node, thin_or_thick):
    matrix = sp.Matrix()
    if thin_or_thick.lower() == "thin":
        
        if node.getType() == "drift teapot":
            
            length = node.getLength()
            matrix = make_drift(length)
            
        elif node.getType() == "quad teapot":
            
            kq = node.getParam("kq")
            length = node.getLength()
            f = 1/(kq*length)
            matrix = make_quad_focus(f)
                
        elif node.getType() == "kick teapot":
            kx = 0
            ky = 0
            length = node.getLength()

            if node.getName() == "IKICKH_A10":
                kx = np.sign(node.getParam("kx"))*kx1
            elif node.getName() == "IKICKH_A11":
                kx = np.sign(node.getParam("kx"))*kx2
            elif node.getName() == "IKICKH_A12":
                kx = np.sign(node.getParam("kx"))*kx3
            elif node.getName() == "IKICKH_A13":
                kx = np.sign(node.getParam("kx"))*kx4
            
            matrix = make_kick_vector(kx)
            
        else:
            print("element not recognized, will add drift instead")
            matrix = Matrix([[1,node.getLength()], [0,1]])
    elif thin_or_thick.lower() == "thick":
        print("I havent implemented this yet")
        matrix = Matrix([[1,0],[0,1]])
            
    return matrix
        
lattice_dict = {}

for n in nodes:
    
    lattice_dict[n.getName()] = make_matrix(n, "thin")
    
def track_through_lattice(nodes, lattice, lattice_dict, kicker_names, end_index, x, xp):
    
    accumulated_matrix = sp.eye(2)
    #vector_kick = sp.eye(2,1)
    xxp = sp.Matrix([x,xp])
    dummy_vector = sp.Matrix()
    dummy_vector2 = sp.Matrix()
    for n in nodes[:end_index+1]:
        pprint(n.getName())
        pprint(accumulated_matrix)
        pprint(lattice_dict[n.getName()])
        if n.getName() in kicker_names:
            
            xxp = accumulated_matrix*xxp
            
            dummy_vector = xxp
            
            accumulated_matrix = sp.eye(2)
            
            kick_vector = lattice_dict[n.getName()]
            
            xxp = xxp + kick_vector
            pprint(xxp)
            
            dummy_vector2 = kick_vector
            
        else:
             accumulated_matrix =  lattice_dict[n.getName()]*accumulated_matrix
    pprint(xxp)
    pprint(accumulated_matrix)
    xxp = accumulated_matrix * xxp
    
    return xxp, accumulated_matrix, dummy_vector, dummy_vector2

           
kicker_names = ["IKICKH_A10", "IKICKH_A11", "IKICKH_A12", "IKICKH_A13"]

xxp, trmatrx, dummy_vector, dummy_vector2 = track_through_lattice(nodes, lattice, lattice_dict, kicker_names, 28, 0, 0)

trMatrix_1 = sp.Matrix([[dummy_vector[0]/kx1, 0],
                        [dummy_vector[1]/kx1, dummy_vector2[1]/kx2]])


trMatrix_total = trmatrx*trMatrix_1
inv_trMatrix_total = trMatrix_total.inv()

print("\ntransfer matrix:\n")
pprint(trMatrix_total)


print("\ninverse matrix:\n")
pprint(inv_trMatrix_total)

x_xpvector = sp.Matrix([x,xp])
kx_kx2_vector = inv_trMatrix_total*x_xpvector
print("\nanalytical solution for kx1, kx2:\n")
pprint(kx_kx2_vector)