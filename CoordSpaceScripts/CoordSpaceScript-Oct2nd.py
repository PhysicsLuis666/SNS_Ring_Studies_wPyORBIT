#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 10:20:54 2025

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
from matplotlib.patches import Polygon, Patch
from matplotlib.path import Path

import time as t
import scipy as sp
from scipy.spatial import ConvexHull
import sympy as sym
import seaborn as sns
sys.path.insert(0,"/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization" )
sys.path.insert(0, "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/modules_to_use")

from orbit.core.bunch import Bunch, BunchTwissAnalysis
import pakagesForOptimizationScripts as pk
np.printoptions(legacy="1.25")

#set up arrays from user directories
#--------------------------------------------------
# image_folder = "/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/SNS_Ring_Studies_wPyORBIT/OptimizingScripts/Sept5thCoImage43"
df_success_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/success_array_4.csv")
df_fail_array = pd.read_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/fail_array_4.csv")
df_combined = pd.concat([df_success_array, df_fail_array], ignore_index=True)

#set up bunch and lattice
#-------------------------------------------------
bunch = Bunch()
pk.do_bunch(bunch, 1.3)
lattice = pk.build_injection_region_lattice()  #teapot ring object 
nodes = lattice.getNodes()
target_node = nodes[28]
Hkicks = pk.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
kicker_strengths = pk.obtain_scaled_kicker_strengths(lattice, "horizontal", Hkicks, [0,0])  
lattice.initialize()
    
#set up array to plot and success condition:
#-----------------------------------------------
x_xp_by_kicker={}
for i, limit in enumerate(kicker_strengths):
    kicker_col=f'k{i+1}'
    abs_limit = abs(limit)
    
    if abs_limit > 0:
        df_ok = df_combined[(abs(df_combined[kicker_col])<= abs_limit)]
    elif abs_limit == 0:
        df_ok = df_combined[df_combined[kicker_col] == 0]
    
    x_xp_by_kicker[kicker_col] = df_ok[["x","xp"]].reset_index(drop=True)


offset_x = [0, 2.5e-4, 6.5e-4, 8.5e-4]
offset_xp = [0, 1e-5, 1e-5, 1e-5]
    
colors = {"k1":"red", "k2":"green", "k3":"blue", "k4":"black"} 

#draw convex hulls
#------------------------------------------------
fig, ax = plt.subplots(figsize=(8,8))   
polygons = {} #store polygons for click detection
legend_patches = []
for j, k in enumerate(["k1", "k2", "k3", "k4"]): 
    df = x_xp_by_kicker[k].copy()
    df["x_offset"] = df["x"] + offset_x[j]
    df["xp_offset"]= df["xp"] + offset_xp[j]
    
    points = df[["x", "xp"]].to_numpy() #this makes the selected columns into numpy arrays
    
    if len(points) >= 3:
         hull = ConvexHull(points) #Convexhull expects 2D numpy arrays as input, other functionality as well
         hull_points = points[hull.vertices] # Indices of points forming the vertices of the convex hull. For 2-D convex hulls, the vertices are in counterclockwise order. For other dimensions, they are in input order.
         hull_points = np.vstack([hull_points, hull_points[0]]) #this closes the shape
        
         poly_patch = Polygon(hull_points, facecolor=colors[k], alpha=.1)
         ax.add_patch(poly_patch)
         # poly = ax.fill(hull_points[:,0],
         #                hull_points[:,1],
         #                color=colors[k],
         #                alpha=.1,
         #                label=k, picker=True)[0]
         ax.plot(hull_points[:,0],
                         hull_points[:,1],
                         color=colors[k],
                         lw=.8)
         
         polygons[k] = Path(hull_points)
         legend_patches.append(Patch(facecolor=colors[k], edgecolor=colors[k], alpha=.3, label=k))
         
    #plt.scatter(df["x_offset"], df["xp_offset"], s=.4, color=colors[k], alpha=.7)
ax.legend(handles=legend_patches)
ax.grid()
ax.set_xlabel("x [m]")
ax.set_ylabel("xp [rad]")
ax.set_title("Coordinate Space Outlines")
plt.tight_layout()

#set up clicker
#--------------------------------------------------
def on_click(event):
    if event.inaxes != ax:
        return
    
    x_click, xp_click = event.xdata, event.ydata
    
    for k, path in polygons.items():
        if path.contains_point([x_click, xp_click]):
            print(f"Click occured in {k}: {x_click, xp_click}")
    
            # Run optimizer & save
            pk.set_optimizer_plot(lattice, "horizontal", target_node, kicker_strengths, x_click, xp_click)
            
            folder_name = "October2CoImageClickImages"
            os.makedirs(folder_name, exist_ok=True)
            fname = f"Closed_Orbit_Plot_x{x_click:.6f}_xp{xp_click:.6f}.png"
            filepath = os.path.join(folder_name, fname)
            
            if os.path.exists(filepath):
                img = mpimg.imread(filepath)
                fig_view, ax_img = plt.subplots(num = "Closed Orbit Viewer", figsize=(6,6))
                ax_img.imshow(img)
                ax_img.set_title(f"Closed Orbit: [x={x_click*1000:.6f} mm, xp={xp_click*1000:.6f} mrad]")
                ax_img.axis("off")
                plt.show(block=False)
                plt.draw()
                plt.pause(0.1)
                
            else:
                print("⚠️ Expected file not found:", filepath)
        
#attach cursor to ax.subplots 
#cursor = mplcursors.cursor(ax, hover=False)
#cursor.connect("add", on_click)

#connect cursor to figure
fig.canvas.mpl_connect("button_press_event", on_click)
plt.show()
