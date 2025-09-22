#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 12:05:47 2025

@author: l5g
"""

#This Script will obtain the transfer matrix anywhere in the lattice 

from contextlib import contextmanager
import os
import time as t
import json
import sys as s 

import math
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import sign
from scipy.optimize import minimize 


from orbit.core import bunch
from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.teapot import DriftTEAPOT, QuadTEAPOT, KickTEAPOT, TEAPOT_Ring, TEAPOT_Lattice
from orbit.injection import TeapotInjectionNode, JohoTransverse, JohoLongitudinal
from orbit.injection import SNSESpreadDist
from orbit.kickernodes import SquareRootWaveform, flatTopWaveform
from orbit.lattice import AccActionsContainer, AccLattice, AccNode, AccNodeBunchTracker
from orbit.diagnostics.TeapotDiagnosticsNode import TeapotBPMSignalNode
from orbit.utils.consts import mass_proton, charge_electron, speed_of_light

import pakagesForOptimizationScripts as pfo
import PlotingFunction as po
import OPTIMIZATIONFunctionPackages as opt

from orbit.core.orbit_utils import Matrix
from orbit.matrix_lattice import BaseMATRIX
from orbit.teapot import TEAPOT_MATRIX_Lattice

from orbit.py_linac.lattice import MarkerLinacNode
from orbit.py_linac.lattice import LinacTrMatricesController

bunch = Bunch()
pfo.do_bunch(bunch, 1.3)
lattice = pfo.build_lattice()
nodes = lattice.getNodes()

drift_a09 = lattice.getNodeForName("DMCV_A09")
drift_b01 = lattice.getNodeForName("DMCV_B01")
#Name:DH_A11B, index:27, type: bend teapot
#Name:DB23, index:28, type: drift teapot
#Name:DDMCV_X01 index:51 type: drift teapot
reordered_lattice = pfo.reorganize_lattice(lattice, drift_b01, drift_a09)
reordered_nodes = reordered_lattice.getNodes()
target_node = reordered_nodes[28] #a drift teapot  'DB23' index 28, at index 27 it is a bend teapot 'DH_A11B'
target_name = target_node.getName()
reordered_lattice.initialize()
horizontal_kickers = pfo.get_kickers(reordered_lattice, "horizontal")
Hkicks = pfo.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
hk_strengths = pfo.obtain_scaled_kicker_strengths(reordered_lattice, "horizontal", Hkicks, [0,0])
pfo.set_kicker_strengths(reordered_lattice, "horizontal", hk_strengths)


for n in reordered_nodes:
    print(f"{n.getName()}, {n.getType()}, index:{reordered_lattice.getNodeIndex(n)}")

def generate_matrix(lattice, nodes, index):
    
    
    return []







































# trMatriciesGenerator = LinacTrMatricesController()
# lattice_trMatrix_nodes = horizontal_kickers
# trMatrixNodes = trMatriciesGenerator.addTrMatrixGenNodes(reordered_lattice, 
#                                                          lattice_trMatrix_nodes,
#                                                          place=MarkerLinacNode.ENTRANCE)

# for trMatrxNode in trMatrixNodes:
#     trMatrxNode.setTwissWeightUse(True,True,True)
    

# reordered_lattice.trackBunch(bunch)


# for trMtrxNode in trMatrixNodes:
#     print ("Transport Matrix Node =",trMtrxNode.getName()," pos[m] = %5.3f"%trMtrxNode.getPosition())
#     trMtrx = trMtrxNode.getTransportMatrix()
#     n = trMtrx.size()[0]
#     m = trMtrx.size()[1]
#     matrix = np.zeros([n,m])
#     for i in range(n):
#         for j in range(m):
#             matrix[i,j] = trMtrx.get(i,j)
#     print(matrix)
#     print ("====================================")
        



