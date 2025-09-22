#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 14:56:05 2025

@author: l5g
"""
#Optimization Package
import random
import sys
import numpy as np
import scipy as sp
from scipy.optimize import minimize
import math as mt
import matplotlib.pyplot as plt
import pakagesForOptimizationScripts as pko

import time as t

from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op
from orbit.core import bunch
from numpy import sign
from orbit.core.bunch import Bunch
from orbit.lattice import AccActionsContainer, AccLattice, AccNode, AccNodeBunchTracker
from orbit.utils.consts import mass_proton, charge_electron, speed_of_light


def optimizing_function(k0, lattice, plane_selection, target_node, phase_space_goal):
    kicker_nodes = pko.get_kickers(lattice, plane_selection)
    k4 = kicker_nodes[-1]
    nodes = lattice.getNodes()
    pko.set_kicker_strengths(lattice, plane_selection, k0)
    
    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(1.300)
    bunch.addParticle(0,0,0,0,0,0)
    
    lattice.trackBunch(bunch)
    
    ps_foil = pko.get_phase_space_at_diagnostic(lattice, nodes, target_node, f"foil_bpm_{target_node.getName()}", "entrance")
    coords_at_kicker4 = pko.get_phase_space_at_diagnostic(lattice, nodes, k4, f"diagnostic_kicks_{k4.getName()}", "exit")
    
    score= 0.0
    for i in range(2):
        score += (ps_foil[i] - phase_space_goal[i])**2 
    score += 1e3*(coords_at_kicker4[0] **2 + coords_at_kicker4[1]**2)
    #print(f"xk4:{coords_at_kicker4[0]}, xpk4:{coords_at_kicker4[1]}")
    
    return score

def get_phase_space_goal(x, xp, lattice, foil_node):
    
    nodes = lattice.getNodes()
    xcenterpos, xcentermom = 0.0486, 0
    ycenterpos, ycentermom = 0.046, 0
   
    pos_foil = pko.determine_distance(nodes, foil_node)

    inj_pt = [pos_foil, xcenterpos]
    inj_ps = [xcenterpos, xcentermom]
    
    phase_space_goal = [inj_ps[0]-x, inj_ps[1] - xp]
    
    return phase_space_goal


def do_optimize(x, xp, lattice, selected_plane, foil_node, kicker_strengths):
    
    phase_space_goal = get_phase_space_goal(x, xp, lattice, foil_node)
    
    res = minimize(optimizing_function, 
                   kicker_strengths, 
                   args=(lattice, selected_plane, foil_node, phase_space_goal),
                   method="Nelder-Mead",
                   options={"xatol":1e-09, "fatol":1e-09, "maxiter":500})
    
    return res.x

    

