#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 12:07:47 2025

@author: l5g

this script will get the csv files from the OptimizationScriptChangeInOutputFile
and then will plot the coord space with color mapings to which kickers fail first

"""

import os 
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time as t
import scipy as sp
import sympy as sym

from orbit.core.bunch import Bunch, BunchTwissAnalysis
from orbit.teapot.teapot_matrix_lattice import TEAPOT_MATRIX_Lattice, TEAPOT_Lattice
from orbit.matrix_lattice import BaseMATRIX, MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccNodeBunchTracker, AccActionsContainer
import pakagesForOptimizationScripts as pk


sys.path.insert(0,"/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization" )
sys.path.insert(0, "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/modules_to_use")


def check_where_kickFail(new_kicker_strengths, old_kicker_strengths):
    
    
    for i, k in enumerate(new_kicker_strengths, start = 1):
        
        limit = old_kicker_strengths[i-1]
        
        if limit > 0:
            if not (0 < k <= limit):
                return i
        elif limit < 0:
            if not (limit <= k < 0):
                return i
        elif limit == 0:
            if k != 0:
                return i
        
    return 0


df_success_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/success_array_2.csv")
df_fail_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/fail_array_2.csv")

bunch = Bunch()
pk.do_bunch(bunch, 1.3)
lattice = pk.build_injection_region_lattice()  #teapot ring object 
Hkicks = pk.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
kicker_strengths = pk.obtain_scaled_kicker_strengths(lattice, "horizontal", Hkicks, [0,0])  

df_fail_array["failed_kickers"] = df_fail_array.apply(lambda row: check_where_kickFail([row.k1, row.k2, row.k3, row.k4], kicker_strengths ), axis=1)

colors = {1:"blue", 2:"yellow", 3:"pink", 4:"orange", 0:"green"}
labels = {1:"Fail kicker 1", 2:"Fail kicker 2", 3:"Fail kicker 3", 4:"Fail kicker 4", 0:"No kicker failed"}
plt.figure(figsize=(6,5))
for k in colors.keys():
    subset = df_fail_array[df_fail_array["failed_kickers"]==k]
    plt.scatter(subset["x"], subset["xp"], s=5,  color=colors[k], alpha=.7, label = labels[k])
    
plt.scatter(df_success_array["x"], df_success_array["xp"], s=10, color="green", alpha=.6, label="kickers success")
plt.xlabel("x [m]")
plt.ylabel("xp [rad]")
plt.grid() 
plt.legend()
plt.title("coordinate kicker space")
plt.show()   
    
    
# plt.scatter(df_fail_array["x"],
#             df_fail_array["xp"],
#             c=df_fail_array["failed_kickers"].map(colors),
#             alpha=.7)
# plt.scatter(df_success_array["x"],
#             df_success_array["xp"],
#             color="green",
#             alpha=.7)
# plt.title("kicker coordinate space")
# plt.xlabel("x [m]")
# plt.ylabel("xp [rad]")

# plt.show()

