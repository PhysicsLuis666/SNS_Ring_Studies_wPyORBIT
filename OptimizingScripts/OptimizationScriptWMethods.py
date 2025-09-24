#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 14:32:12 2025

@author: l5g
"""
import os
import time as t
import json
import sys as s 
from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op
import math
import random
import sys
import pandas as pd
from orbit.core import bunch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from numpy import sign
from scipy.optimize import minimize 
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

e = charge_electron
m_p = mass_proton
c = speed_of_light
bunch = Bunch()
bunch.mass(m_p)
bunch.charge(e)
bunch.getSyncParticle().kinEnergy(1.3)
bunch.addParticle(0,0,0,0,0,0)
lattice = pfo.build_lattice()
nodes = lattice.getNodes()
print(nodes[0].getName())

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

hkickers = pfo.get_kickers(reordered_lattice, "horizontal")
for node in reordered_nodes:
    if node.getName() == "DMCV_A09":
        pfo.place_diagnostic_node(node, f"diagnostic_{node.getName()}", "exit")
    elif node.getName() == "DMCV_B01":
        pfo.place_diagnostic_node(node, f"diagnostic_{node.getName()}", "entrance")
    elif node == hkickers[0]:
        pfo.place_diagnostic_node(node, f"diagnostic_kicks_{node.getName()}", "entrance")
    elif node in [hkickers[1],hkickers[2]]:
        pfo.place_diagnostic_node(node, f"diagnostic_kicks_{node.getName()}", "entrance")
    elif node == hkickers[3]:
        pfo.place_diagnostic_node(node, f"diagnostic_kicks_{node.getName()}", "exit")
    elif node == target_node:
        pfo.place_diagnostic_node(target_node, f"foil_bpm_{target_name}", "entrance")
    else:
        pfo.place_diagnostic_node(node, f"diagnostic_{node.getName()}", "entrance")

Hkicks = pfo.max_angle_deflection(1.3, 1, 1400, 1600, "horizontal")
kicker_strengths = pfo.obtain_scaled_kicker_strengths(reordered_lattice, "horizontal", Hkicks, [0,0])    
pfo.set_kicker_strengths(reordered_lattice, "horizontal", kicker_strengths)

k_success = []
x_coord_success = []
xp_coord_success = []
x_coord_fail = []
xp_coord_fail = []
k_fail = []
sx_coord_foil = []
sxp_coord_foil = []
sxk4_coord = []
sxpk4_coord = []
fx_coord_foil = []
fxp_coord_foil = []
fxk4_coord = []
fxpk4_coord = []
sphase_Goal = []
fphase_Goal = []
for x in np.linspace(0.006, 0.0486, 20):   # 50 points
    for xp in np.linspace(-0.0009, 0.0009, 20):  # 50 points
        new_kicker_strengths = opt.do_optimize(x, xp, reordered_lattice, "horizontal", target_node, kicker_strengths)
        success = pfo.check_lims(kicker_strengths, new_kicker_strengths)
        coord_foil = pfo.get_phaseSpaceCoords(reordered_lattice, new_kicker_strengths, target_node,  f"foil_bpm_{target_node.getName()}", "entrance")
        coord_kicker4 = pfo.get_phaseSpaceCoords(reordered_lattice, new_kicker_strengths, hkickers[3], f"diagnostic_kicks_{hkickers[3].getName()}", "exit")
        phaseSpace_goal = opt.get_phase_space_goal(x, xp, reordered_lattice, target_node)

        if success == True:
            print("-------success-------")
            print(f"[hase_space_goal:{phaseSpace_goal}")
            print(f"x:{x},xp:{xp}, k:{new_kicker_strengths}")
            print(f"xk4:{coord_kicker4[0]}, xpk4:{coord_kicker4[1]}")
            print(f"coords at foil: x:{coord_foil[0]}, xp:{coord_foil[1]}")
            po.plot_w_new_strengths(reordered_lattice, hkickers, new_kicker_strengths, "Sept5thCoImages")
            k_success.append(new_kicker_strengths)
            x_coord_success.append(x)
            xp_coord_success.append(xp)
            sx_coord_foil.append(coord_foil[0])
            sxp_coord_foil.append(coord_foil[1])
            sxk4_coord.append(coord_kicker4[0])
            sxpk4_coord.append(coord_kicker4[1])
            sphase_Goal.append(phaseSpace_goal)
            
        elif success == False:
            print("-------fail--------")
            print(f"[hase_space_goal:{phaseSpace_goal}")
            print(f"x:{x},xp:{xp}, k:{new_kicker_strengths}")
            print(f"xk4:{coord_kicker4[0]}, xpk4:{coord_kicker4[1]}")
            print(f"coords at foil: x:{coord_foil[0]}, xp:{coord_foil[1]}")
            k_fail.append(new_kicker_strengths)
            x_coord_fail.append(x)
            xp_coord_fail.append(xp)
            fx_coord_foil.append(coord_foil[0])
            fxp_coord_foil.append(coord_foil[1])
            fxk4_coord.append(coord_kicker4[0])
            fxpk4_coord.append(coord_kicker4[1])
            fphase_Goal.append(phaseSpace_goal)
        
print(f"successfull runs: {len(k_success)}")
print(f"failed runs: {len(k_fail)}")
plt.scatter(x_coord_success, xp_coord_success, color="green")
plt.scatter(x_coord_fail, xp_coord_fail, color="red")
plt.xlabel("x values [m]")
plt.ylabel("xp values[rad]")
#plt.xlim(.0099, .0486)
os.makedirs("x-xpSpace", exist_ok=True)  # Make sure directory exists
output_path = os.path.join("x-xpSpace", "allowedx-x'values.png")
plt.savefig(output_path, dpi=300, bbox_inches="tight")

#k_success = np.array(k_success)
success_array_coords = np.stack([x_coord_success, xp_coord_success, sx_coord_foil, sxp_coord_foil, sxk4_coord, sxpk4_coord], axis=1)
fail_array_coords = np.stack([x_coord_fail, xp_coord_fail, fx_coord_foil, fxp_coord_foil, fxk4_coord, fxpk4_coord], axis=1)

cols = ["x", "xp", "x_foil", "xp_foil", "x_kicker4", "xp_kicker4"]

df_success = pd.DataFrame(success_array_coords, columns=cols)
df_fail    = pd.DataFrame(fail_array_coords, columns=cols)

df_success.to_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/success_coords.csv", index=False)
df_fail.to_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/fail_coords.csv", index=False)

dgoal_success = pd.DataFrame(sphase_Goal, columns=["x","xp"])
dgoal_fail = pd.DataFrame(fphase_Goal, columns=["x","xp"])

dgoal_success.to_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/success_ps_goal.csv", index=False)
dgoal_fail.to_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/fail_ps_goal.csv", index=False)

dk_success = pd.DataFrame(k_success, columns=["k1","k2","k3","k4"])
dk_fails = pd.DataFrame(k_fail, columns=["k1","k2","k3","k4"])

dk_success.to_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/k_success.csv", index=False)
dk_fails.to_csv("/Users/l5g/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/sns-ring-optimizationStudies-2025/August2025OptimizationScripts/September2025OptimizationScripts/DataOptimization/k_fails.csv", index=False)