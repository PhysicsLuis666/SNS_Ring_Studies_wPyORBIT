#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 15:12:20 2025

@author: l5g
"""
import os 
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time as t
import scipy as sp
from scipy.spatial import ConvexHull
import sympy as sym
import seaborn as sns
sys.path.insert(0,"/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization" )
sys.path.insert(0, "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/modules_to_use")

from orbit.core.bunch import Bunch, BunchTwissAnalysis
from orbit.teapot.teapot_matrix_lattice import TEAPOT_MATRIX_Lattice, TEAPOT_Lattice
from orbit.matrix_lattice import BaseMATRIX, MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccNodeBunchTracker, AccActionsContainer
import pakagesForOptimizationScripts as pk





df_success_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/success_array_3.csv")
df_fail_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/fail_array_3.csv")
df_combined = pd.concat([df_success_array, df_fail_array], ignore_index=True)


bunch = Bunch()
pk.do_bunch(bunch, 1.3)
lattice = pk.build_injection_region_lattice()  #teapot ring object 
Hkicks = pk.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
kicker_strengths = pk.obtain_scaled_kicker_strengths(lattice, "horizontal", Hkicks, [0,0])  

x_xp_by_kicker={}

for i, limit in enumerate(kicker_strengths):
    kicker_col=f'k{i+1}'
    
    if limit > 0:
        df_ok = df_combined[(df_combined[kicker_col] > 0) & (df_combined[kicker_col] <= limit)]
    elif limit < 0:
        df_ok = df_combined[(df_combined[kicker_col] >= limit) & (df_combined[kicker_col] < 0)]
    else:
        df_ok = df_combined[df_combined[kicker_col] == 0]
    
    x_xp_by_kicker[kicker_col] = df_ok[["x","xp"]].reset_index(drop=True)

#print(x_xp_by_kicker["k1"])
#print(x_xp_by_kicker["k2"])

offset_x = [0, 2.5e-4, 6.5e-4, 8.5e-4]
offset_xp = [0, 1e-5, 1e-5, 1e-5]
# for j, kicker in enumerate(['k1','k2','k3','k4']):
#    x_xp_by_kicker[kicker]["x_offset"] = x_xp_by_kicker[kicker]["x"] + offset_x[j]
#    x_xp_by_kicker[kicker]["xp_offset"] = x_xp_by_kicker[kicker]["xp"] + offset_xp[j]
    
    
colors = {"k1":"red", "k2":"green", "k3":"blue", "k4":"black"} 

plt.figure(figsize=(8,8))   

for j, k in enumerate(["k1", "k2", "k3", "k4"]): 
    df = x_xp_by_kicker[k].copy()
    df["x_offset"] = df["x"] + offset_x[j]
    df["xp_offset"]= df["xp"] + offset_xp[j]
    
    points = df[["x_offset", "xp_offset"]].to_numpy() #this makes the selected columns into numpy arrays
    
    if len(points) >= 3:
        hull = ConvexHull(points) #Convexhull expects 2D numpy arrays as input, other functionality as well
        hull_points = points[hull.vertices] # Indices of points forming the vertices of the convex hull. For 2-D convex hulls, the vertices are in counterclockwise order. For other dimensions, they are in input order.
        hull_points = np.vstack([hull_points, hull_points[0]]) #this closes the shape
        
        plt.fill(hull_points[:,0], hull_points[:,1], color=colors[k], alpha=.1, label=k)
        plt.plot(hull_points[:,0], hull_points[:,1], color=colors[k], lw=1)
    

plt.legend()
plt.grid()
plt.xlabel("x [m]")
plt.ylabel("xp [rad]")
plt.title("Coordinate Space Outlines")
plt.xlim(.0055,.05)
plt.ylim(-.001, .001)
plt.show()
    
    
    

# plt.scatter(x_xp_by_kicker["k1"]["x_offset"], x_xp_by_kicker["k1"]["xp_offset"], s=5, color="red", alpha=.4, label="k1")
# plt.scatter(x_xp_by_kicker["k2"]["x_offset"], x_xp_by_kicker["k2"]["xp_offset"], s=4, color="green", alpha=.7, label="k2")
# plt.scatter(x_xp_by_kicker["k3"]["x_offset"], x_xp_by_kicker["k3"]["xp_offset"], s=3, color="blue", alpha=.7, label="k3")
# plt.scatter(x_xp_by_kicker["k4"]["x_offset"], x_xp_by_kicker["k4"]["xp_offset"], s=3, color="black", alpha=.7, label="k4")
# plt.xlabel("x [m]")
# plt.xlim(.0055, .05)
# plt.ylim(-.001, .001)
# plt.ylabel("xp [rad]")
# plt.title("Coordinate Space")
# plt.legend()
# plt.grid()
# plt.show()

#sns.kdeplot(x=x_xp_by_kicker['k1']['x'],
           # y=x_xp_by_kicker['k2']['xp'],
           # levels=5,   fills=False, color="red", linewidth=2)




