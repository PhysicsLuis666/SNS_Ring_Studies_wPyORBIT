#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 14:37:10 2025

@author: l5g

This script will obtain the transport matrix of a lattice using andries example
LinacTRMatrix method in PyORBIT 3

This script will not work with the Ring TEAPOT Class, personally asked
andrie and he cleard it up 

"""

import os 
import sys
import time as t
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import sympy as syp
sys.path.insert(0, "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/modules_to_use")
from orbit.core.bunch import Bunch, BunchTwissAnalysis
from orbit.teapot.teapot_matrix_lattice import TEAPOT_MATRIX_Lattice, TEAPOT_Lattice
from orbit.matrix_lattice import BaseMATRIX, MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccNodeBunchTracker, AccActionsContainer
from orbit.py_linac.lattice import MarkerLinacNode, LinacTrMatricesController
import pakagesForOptimizationScripts as pk

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

trMatrixGenerator = LinacTrMatricesController() 
lattice_trMatrix_Nodes = hkicker_nodes
trMatrixNodes = trMatrixGenerator.addTrMatrixGenNodes(lattice, lattice_trMatrix_Nodes, place=MarkerLinacNode.EXIT)
print(trMatrixNodes)

# for mtrx in trMatrixNodes:
#     mtrx.setTwissWeightUse(True,True,True)
    
bunch = Bunch()
    
bunch1.copyBunchTo(bunch)


lattice.trackBunch(bunch)

for trmx in trMatrixNodes:
    print(f"transport matrix node:{trmx.getName()}, at pos[m] = {trmx.getPosition()} ")
    trMatrix = trmx.getTransportMatrix()
    tr_matrix = pk.returnMtrx(trMatrix)
    print(tr_matrix)
    print ("====================================")


