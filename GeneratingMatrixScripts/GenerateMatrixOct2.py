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
ky1,ky2,ky3,ky4 = sp.symbols(r"ky1,ky2,ky3,ky4")
f1,f2,f3,f4 = sp.symbols(r"f1,f2,f3,f4")
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10 = sp.symbols(r"l1, l2, l3, l4, l5, l6, l7, l8, l9, l10")

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
    
    # if k != 0:
    #     vector = sp.Matrix([0, k])
    # else:
    #     vector = sp.eye(2)
    return sp.Matrix([0,k])
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
pk.set_kicker_strengths(lattice, "horizontal", hk_strengths)
lattice.initialize()   
def make_matrix(node, thin_or_thick):
    matrix = sp.Matrix()
    if thin_or_thick.lower() == "thin":
        
        if node.getType() == "drift teapot":
            
            length = node.getLength()
            matrix = make_drift(length)
            matrix = [matrix]
            
        elif node.getType() == "quad teapot":
            
            kq = node.getParam("kq")
            length = node.getLength()
            f = 1/(kq*length)
            drift_matrix = make_drift(length/2)
            quad_matrix = make_quad_focus(f)
            matrix = [drift_matrix, quad_matrix, drift_matrix]
                
        elif node.getType() == "kick teapot":
            kx = 0
            ky = 0
            
            length = node.getLength() 
                
            if length != 0:
                l = length/2
                drift_matrix_kicks = make_drift(l)
            
            if node.getName() == "IKICKH_A10":
                print(f"here is the kicker parameters") 
                print(f"k1:{n.getParam('kx')}")
                if np.sign(node.getParam("kx")) < 0:
                    kx = -kx1
                else:
                    kx = kx1
               
            elif node.getName() == "IKICKH_A11":
                print(f"here is the kicker parameters")
                print(f"k2:{n.getParam('kx')}")
                if np.sign(node.getParam("kx")) < 0:
                    kx = -kx2
                else:
                    kx = kx2
               
            elif node.getName() == "IKICKH_A12":
                print(f"here is the kicker parameters")
                print(f"k3:{n.getParam('kx')}")
                if np.sign(node.getParam("kx")) < 0:
                    kx = -kx3
                else:
                    kx = kx3 
            elif node.getName() == "IKICKH_A13":
                print(f"here is the kicker parameters")
                print(f"k4:{n.getParam('kx')}")
                if np.sign(node.getParam("kx")) < 0:
                    kx = -kx4
                else:
                    kx = kx4
          
                 
            kicker_vector = make_kick_vector(kx)
            
        
            if length != 0  n.getName() in ['IKICKV_A10', 'IKICKV_A11', 'IKICKV_A12', 'IKICKV_A13']:
                    matrix = [drift_matrix_kicks, make_kick_vector(0), drift_matrix_kicks]
            elif length != 0:
                    matrix = [drift_matrix_kicks, kicker_vector, drift_matrix_kicks]
            else:
                matrix = make_kick_vector(kx)
                matrix = [matrix]
             

            
        else:
            print("element not recognized, will add drift instead")
            matrix = Matrix([[1,node.getLength()], [0,1]])
            matrix = [matrix]
    elif thin_or_thick.lower() == "thick":
        print("I havent implemented this yet")
        matrix = Matrix([[1,0],[0,1]])
        matrix = [matrix]
            
    return matrix
        
lattice_dict = {}

for n in nodes:
    lattice_dict[n.getName()] = make_matrix(n, "thin")
    
def track_through_lattice(nodes, lattice, lattice_dict,start_index, end_index, x, xp):
    
    accumulated_matrix = sp.eye(2)
    #vector_kick = sp.eye(2,1)
    xxp = sp.Matrix([x,xp])
    kicks_applied = []
    states = []
    
    for n in nodes[start_index:end_index+1]:
        for i in lattice_dict[n.getName()]:
            print("currently on")
            print(f"node:{n.getName()} with matrix:")
            pprint(i)
            if i.shape == (2,2):
                    accumulated_matrix = i*accumulated_matrix
                    
            elif i.shape == (2,1):
                xxp = accumulated_matrix*xxp
                
                accumulated_matrix = sp.eye(2)
                
                xxp = xxp + i
                
                kicks_applied.append(i)
                states.append(xxp) #stores the kick vector we currently are on
            print("accumulated matrix:")
            pprint(accumulated_matrix)
            print("current coordinate xxp:")
            pprint(xxp)
        
    #pprint(xxp)
    #pprint(accumulated_matrix)
    
    xxp = accumulated_matrix * xxp
    
    return xxp, accumulated_matrix, kicks_applied, states
    

print("starting first half")
xxp, trmatrixd5tod3, kicks, states = track_through_lattice(nodes, lattice, lattice_dict, 0, 28, 0, 0)



#second_part
print("starting second half")
xxp2, accumlated_matrix2, kicks_applied, states2 = track_through_lattice(nodes, lattice, lattice_dict, 28, 53, .014669, -.000305)   



