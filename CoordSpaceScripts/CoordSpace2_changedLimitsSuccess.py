#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 14:56:07 2025

@author: l5g
"""

import os 
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

import mplcursors
import matplotlib.image as mpimg
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

np.printoptions(legacy="1.25")

# image_folder = "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/OptimizingScripts/Sept5thCoImage43"
df_success_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/success_array_4.csv")
df_fail_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/fail_array_4.csv")
df_combined = pd.concat([df_success_array, df_fail_array], ignore_index=True)

# #Regex: match [x,xp] with comma seperator
# patern = re.compile(r"Closed_Orbit_Plot_\[([-\d.eE]+),\s*([-\d.eE]+)\]\.png")
# coord_to_file = {}

# for fname in os.listdir(image_folder):
#     m = patern.match(fname)
    
#     if m:
#         x, xp = float(m.group(1)), float(m.group(2))
#         coord_to_file[(x,xp)] = os.path.join(image_folder, fname)

# def find_closest_file(x,xp):
#     coords = np.array(list(coord_to_file.keys()))
#     idx = np.argmin(np.linalg.norm(coords - np.array([x,xp]), axis=1))
#     return list(coord_to_file.values())[idx]


bunch = Bunch()
pk.do_bunch(bunch, 1.3)
lattice = pk.build_injection_region_lattice()  #teapot ring object 
nodes = lattice.getNodes()
target_node = nodes[28]
Hkicks = pk.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
kicker_strengths = pk.obtain_scaled_kicker_strengths(lattice, "horizontal", Hkicks, [0,0])  

x_xp_by_kicker={}
for i, limit in enumerate(kicker_strengths):
    kicker_col=f'k{i+1}'
    abs_limit = abs(limit)
    
    if abs_limit > 0:
        df_ok = df_combined[(abs(df_combined[kicker_col])<= abs_limit)]
    elif abs_limit == 0:
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

fig, ax = plt.subplots(figsize=(8,8))   

for j, k in enumerate(["k1", "k2", "k3", "k4"]): 
    df = x_xp_by_kicker[k].copy()
    df["x_offset"] = df["x"] + offset_x[j]
    df["xp_offset"]= df["xp"] + offset_xp[j]
    
    points = df[["x", "xp"]].to_numpy() #this makes the selected columns into numpy arrays
    
    if len(points) >= 3:
         hull = ConvexHull(points) #Convexhull expects 2D numpy arrays as input, other functionality as well
         hull_points = points[hull.vertices] # Indices of points forming the vertices of the convex hull. For 2-D convex hulls, the vertices are in counterclockwise order. For other dimensions, they are in input order.
         hull_points = np.vstack([hull_points, hull_points[0]]) #this closes the shape
        
         ax.fill(hull_points[:,0], hull_points[:,1], color=colors[k], alpha=.1, label=k)
         ax.plot(hull_points[:,0], hull_points[:,1], color=colors[k], lw=.8)
    
    #plt.scatter(df["x_offset"], df["xp_offset"], s=.4, color=colors[k], alpha=.7)

ax.legend()
ax.grid()
ax.set_xlabel("x [m]")
ax.set_ylabel("xp [rad]")
ax.set_title("Coordinate Space Outlines")
#plt.xlim(.0055,.05)
#plt.ylim(-.001, .001)
plt.show()

# def on_click(sel):
#     x_click, xp_click = sel.target
    
#     fname = find_closest_file(x_click, xp_click)
    
#     #open image in new figure
#     img = mpimg.imread(fname)
#     plt.figure()
#     plt.imshow(img)
#     plt.title(f"Closest image to [{x_click*1000,xp_click*1000}] ")
#     plt.axis("off")
#     plt.show()
    
# #activate cursor

# cursor = mplcursors.cursor(ax, hover=False)
# cursor.connect("add", on_click)
# plt.show()