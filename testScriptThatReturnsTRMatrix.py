#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 18:01:10 2025

@author: l5g
"""
from sympy import init_printing, pprint
init_printing(use_unicode=True)
from sympy import symbols
from sympy.matrices import Matrix
from sympy import eye, simplify, MatrixSymbol
import os
import time
import numpy as np
import matplotlib.pyplot as plt

from orbit.teapot import DriftTEAPOT, QuadTEAPOT, KickTEAPOT

from orbit.core import bunch
from orbit.lattice import AccLattice, AccActionsContainer, AccNode
from orbit.core.bunch import Bunch, BunchTuneAnalysis, BunchTwissAnalysis
from orbit.teapot import TEAPOT_Lattice

from orbit.utils.consts import mass_proton, charge_electron, speed_of_light

lattice = TEAPOT_Lattice()

drift1 = DriftTEAPOT("drift1")
drift1.setLength(1)
quad1 = QuadTEAPOT("test quad1")
quad1.addParam("kq", 1.4)
quad1.setLength(1)
kick_node = KickTEAPOT("kicker1")
kick_node.setParam("kx", 3.0e-3)
lattice.addNode(drift1)
lattice.addNode(quad1)
lattice.addNode(kick_node)

lattice.initialize()
nodes = lattice.getNodes()
print(nodes)


def generate_matrix(node):
    
    matrix = Matrix([])
    
    if node.getType() == "drift teapot":
        matrix = Matrix([[1, nodes[0].getLength()], [0,1]])
    elif node.getType() == "quad teapot":
        matrix = Matrix([[1, 0],[-1/(nodes[1].getParam("kq")*nodes[1].getLength()),1]])
    elif node.getType() == "kick teapot":
        matrix = Matrix([[0],[node.getParam("kx")]])
        
    return matrix

def generate_transferMatrix(nodes):
    
    matrix = generate_matrix(nodes[0])
    vector = Matrix([])
    for n in nodes[1:]:
        pprint(f"currently in {n}")
        if n.getType() == "kick teapot":
            vector += generate_matrix(n)
        else:
            matrix *= generate_matrix(n)
        print(matrix)
        
    
    return matrix
    
reversed_nodes = list(reversed(nodes))

print(reversed_nodes)

trMatrix = generate_transferMatrix(reversed_nodes)

pprint(trMatrix)

# drift_matrix = generate_matrix(nodes[0])
# quad_matrix = generate_matrix(nodes[1])
# kick_vector = generate_matrix(nodes[2])
# pprint(kick_vector)
















# def generate_matrix(node):
#     matrix = np.zeros([2,2])
    
#     if node.getType() == "drift teapot":
#         matrix[0,0] = 1
#         matrix[0,1] = node.getLength()
#         matrix[1,0] = 0
#         matrix[1,1] = 1
#     elif node.getType() == "quad teapot":
#         matrix[0,0] = 1
#         matrix[0,1] = 0
#         matrix[1,0] = node.getParam("kq")*node.getLength()
#         matrix[1,1] = 1
#     else:
#         print("I haven't added any other type of matrix")
        
#     return matrix
    
# drift_mtr = generate_matrix(nodes[0])
