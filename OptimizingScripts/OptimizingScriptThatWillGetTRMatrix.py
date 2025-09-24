#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 17:22:46 2025

@author: l5g
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

import pakagesForOptimizationScripts as pk

bunch = Bunch()
pk.do_bunch(bunch, 1.3)
lattice = pk.build_injection_region_lattice()  #teapot ring object 




