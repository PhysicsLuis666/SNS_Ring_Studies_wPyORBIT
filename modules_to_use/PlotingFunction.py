#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 11:45:04 2025

@author: l5g
"""

import matplotlib.pyplot as plt
import numpy as np
import pakagesForOptimizationScripts as pfo
import OPTIMIZATIONFunctionPackages as OP
from scipy.optimize import minimize
from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op
from orbit.core import bunch
from numpy import sign
from orbit.core.bunch import Bunch
from orbit.lattice import AccActionsContainer, AccLattice, AccNode, AccNodeBunchTracker
from orbit.utils.consts import mass_proton, charge_electron, speed_of_light


import os
import matplotlib.pyplot as plt

def plot(lattice, kickers, output_dir):
    nodes = lattice.getNodes()
    x, xp, y, yp, z = pfo.get_coords(lattice)

    x_k1, xp_k1, y_k1, yp_k1, z_k1 = pfo.get_phase_space_at_diagnostic(
        lattice, nodes, kickers[0], "diagnostic_kicks_IKICKH_A10", "entrance")
    x_k2, xp_k2, y_k2, yp_k2, z_k2 = pfo.get_phase_space_at_diagnostic(
        lattice, nodes, kickers[1], "diagnostic_kicks_IKICKH_A11", "entrance")
    x_k3, xp_k3, y_k3, yp_k3, z_k3 = pfo.get_phase_space_at_diagnostic(
        lattice, nodes, kickers[2], "diagnostic_kicks_IKICKH_A12", "entrance")
    x_k4, xp_k4, y_k4, yp_k4, z_k4 = pfo.get_phase_space_at_diagnostic(
        lattice, nodes, kickers[3], "diagnostic_kicks_IKICKH_A13", "exit")

    # Convert x to mm
    x1 = [i * 1000 for i in x]

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(z, x1, color="red")
    plt.plot(z_k1, x_k1 * 1000, "|", markersize=9, color="blue")
    plt.plot(z_k2, x_k2 * 1000, "|", markersize=9, color="blue")
    plt.plot(z_k3, x_k3 * 1000, "|", markersize=9, color="blue")
    plt.plot(z_k4, x_k4 * 1000, "|", markersize=9, color="blue",
             label=f"k4 x:{x_k4}, xp:{xp_k4}, z:{z_k4}")

    plt.grid()
    plt.legend(fontsize=6)

    # --- Save output ---
    os.makedirs(output_dir, exist_ok=True)  # Make sure directory exists
    output_path = os.path.join(output_dir, "Closed_Orbit_plot.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")

    plt.close()  # close instead of show to avoid popup
    print(f"Plot saved to: {output_path}")

def plot_w_new_strengths(lattice, kickers, kicker_strengths, output_dir):
    
    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(1.300)
    bunch.addParticle(0,0,0,0,0,0)
    
    nodes = lattice.getNodes()
    target_node = nodes[28]
    pfo.set_kicker_strengths(lattice, "horizontal", kicker_strengths)
    lattice.trackBunch(bunch)
    
    x,xp,y,yp,z = pfo.get_coords(lattice)
    x_k1, xp_k1, y_k1, yp_k1, z_k1 = pfo.get_phase_space_at_diagnostic(lattice, nodes, kickers[0], 
                                                                       "diagnostic_kicks_IKICKH_A10",
                                                                       "entrance")
    x_k2, xp_k2, y_k2, yp_k2, z_k2 = pfo.get_phase_space_at_diagnostic(lattice, nodes, kickers[1],
                                                                       "diagnostic_kicks_IKICKH_A11", 
                                                                       "entrance")
    x_k3, xp_k3, y_k3, yp_k3, z_k3 = pfo.get_phase_space_at_diagnostic(lattice, nodes, kickers[2],
                                                                       "diagnostic_kicks_IKICKH_A12", 
                                                                       "entrance")
    x_k4, xp_k4, y_k4, yp_k4, z_k4 = pfo.get_phase_space_at_diagnostic(lattice, nodes, kickers[3],
                                                                       "diagnostic_kicks_IKICKH_A13", 
                                                                       "exit")
    x_f, xp_f, y_f, yp_f, z_f = pfo.get_phaseSpaceCoords(lattice, kicker_strengths, target_node,  f"foil_bpm_{target_node.getName()}", "entrance")
    
    #x_b01, xp_b01, y_b01, yp_b01, z_b01 = pfo.get_phase_space_at_diagnostic(lattice, nodes, drift_b01, "diagnostic_DMCV_B01", "entrance")
    x1 = [i*1000 for i in x]
    plt.plot(z, x1, color="red")
    plt.plot(z_k1, x_k1*1000, "|", markersize=9, color="blue")
    plt.plot(z_k2, x_k2*1000, "|", markersize=9, color="blue")
    plt.plot(z_k3, x_k3*1000, "|", markersize=9, color="blue")
    plt.plot(z_k4, x_k4*1000, "|", markersize=9, color="blue", label= f"k4 x:{x_k4}, xp:{xp_k4}, k:{kicker_strengths}")
    plt.plot(z_f, x_f*1000, "o", markersize=6, color="red", label=f"diagnostic on CO RF : x:{x_f*1000}, xp:{xp_f*1000}")
    plt.plot(z_f, .0486*1000, "s", markersize=6, color="orange", label=f"Foil Location x:{.0486*1000}, pos:{z_f}")
    plt.xlim(0,lattice.getLength())
    plt.ylim(-1, 55)
    plt.xlabel("pos [m]")
    plt.ylabel("x [mm]")
    #plt.plot(z_b01, x_b01*1000, "s", markersize=4, color="purple", label=f"B01 x:{x_b01}, xp:{xp_b01}, z:{z_b01}")
    plt.grid()
    plt.legend(fontsize=6)
     
    # --- Save output ---
    os.makedirs(output_dir, exist_ok=True)  # Make sure directory exists
    output_path = os.path.join(output_dir, f"Closed_Orbit_Plot_{kicker_strengths}.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
     
    plt.close()  # close instead of show to avoid popup
    #print(f"Plot saved to: {output_path}")
    

    
    
    