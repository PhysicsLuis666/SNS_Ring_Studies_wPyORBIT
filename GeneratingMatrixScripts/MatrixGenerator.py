#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 14:45:16 2025

@author: l5g
"""
import os 
import sys
import time as t
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from sympy import init_printing, pprint
init_printing(use_unicode=True)
from sympy import symbols as s
from sympy.matrices import Matrix
from sympy import eye, simplify, MatrixSymbol
import sympy as syp

sys.path.insert(0, "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/modules_to_use")
from orbit.core.bunch import Bunch, BunchTwissAnalysis
from orbit.teapot.teapot_matrix_lattice import TEAPOT_MATRIX_Lattice, TEAPOT_Lattice
from orbit.matrix_lattice import BaseMATRIX, MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccNodeBunchTracker, AccActionsContainer
from orbit.py_linac.lattice import MarkerLinacNode, LinacTrMatricesController
import pakagesForOptimizationScripts as pk

kx, ky, f, l, 
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

for n in nodes:
    print(f"{n.getName()}:{n.getType()}, length:{n.getLength()}")
    if n.getType() == "quad teapot":
        print("in quad")
        print(f"k:{n.getParam('kq')}, length:{n.getLength()}")
    if n.getType() == "kick teapot":
        print("in kick")
        print(f"kx:{n.getParam('kx')}")
        print(f"ky:{n.getParam('ky')}")
        print(f"length:{n.getLength()}")
    if n.getType() == "drift teapot":
        print("in drift")

def make_matrix(node, thin_or_thick):
    
    matrix = Matrix()
    length = 0
    if thin_or_thick.lower() == "thin":
        
        if node.getType() == "drift teapot":
            
            length = node.getLength()
            matrix = Matrix([[1,length],[0,1]])
            
        elif node.getType() == "quad teapot":
            
            kq = node.getParam("kq")
            length = node.getLength()
            f = 1/kq*length
            matrix = Matrix([[1,0], [-1/f, 1]])
                
        elif n.getType() == "kick teapot":
            kx = 0
            ky = 0

            if node.getName() in ["IKICKH_A10" , "IKICKH_A11" , "IKICKH_A12" ,"IKICKH_A13"]:
                kx = node.getParam("kx")
            
                matrix = Matrix([[0],[kx]])
            
        else:
            print("element not recognized")
            matrix = Matrix([[1,0], [0,1]])
    elif thin_or_thick.lower() == "thick":
        print("I havent implemented this yet")
        matrix = Matrix([[1,0],[0,1]])
            
        
    return matrix
            
        
lattice_dict = {}

for n in nodes:
    
    lattice_dict[n.getName()] = make_matrix(n, "thin")


# import numpy as np

# I2 = np.eye(2)

# def make_element_matrix(node):
#     """Return (M, b) for the element node (2x2, 2x1)."""
#     t = node.getType()
#     L = node.getLength()
#     if t == "drift teapot":
#         M = np.array([[1.0, L],
#                       [0.0, 1.0]])
#         b = np.zeros((2,1))
#     elif t == "quad teapot":
#         kq = node.getParam("kq")
#         # thin-lens approx: focal strength = kq * L, 1/f = kq*L
#         if abs(kq*L) < 1e-16:
#             M = np.array([[1.0, L],
#                           [0.0, 1.0]])
#             b = np.zeros((2,1))
#         else:
#             f = 1.0/(kq*L)
#             M = np.array([[1.0, 0.0],
#                           [-1.0/f, 1.0]])
#             b = np.zeros((2,1))
#     elif t == "kick teapot":
#         kx = node.getParam("kx") if node.getParam("kx") is not None else 0.0
#         M = np.eye(2)
#         b = np.array([[0.0],[kx]])   # bias: x unchanged, x' += kx
#     else:
#         M = np.eye(2)
#         b = np.zeros((2,1))
#     return M, b

# def compose_sequence(elements):
#     """elements is list of (M,b) in order from entrance->exit.
#        Returns (M_tot, b_tot).
#     """
#     M_tot = np.eye(2)
#     b_tot = np.zeros((2,1))
#     # We must apply elements in order: v -> M1 v + b1 -> M2(...)
#     for M, b in elements:
#         # new total:
#         b_tot = M @ b_tot + b
#         M_tot = M @ M_tot
#     return M_tot, b_tot

# # example usage:
# elements = []
# for node in nodes:           # your lattice node list
#     M, b = make_element_matrix(node)
#     elements.append((M,b))

# Mtot, btot = compose_sequence(elements)
# # To propagate an input vector v0:
# v0 = np.array([[0],[0]])
# vfinal = Mtot @ v0 + btot

        
