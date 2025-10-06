#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 11:46:15 2025

@author: l5g
"""

#import packages
import os
import time as t
from contextlib import contextmanager
import json
import sys as s 
from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op
import math
import random
import sys
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
import OPTIMIZATIONFunctionPackages as ptf
import PlotingFunction as plo
#-----------------------------------------------------------------------------

def build_lattice():
    lattice = TEAPOT_Ring()
   
    lattice.readMAD("/Users/l5g/sns_ring_model/scripts/full_injection_benchmark/inputs/sns_ring_mad.lattice", "RINGINJ")
   
    return lattice 

def get_node_index(lattice, node):
    index = lattice.getNodeIndex(node)
    return index


def reorganize_lattice(lattice, start_node, end_node):
    
    nodes = lattice.getNodes()
    index_start = get_node_index(lattice, start_node)
    index_end = get_node_index(lattice, end_node)
                                                        #index B01 24
    back_end_nodes = nodes[0:index_start+1] #should go from 0 to B01
    front_end_nodes = nodes[index_end:] #this goes from A09 to end 
    
    reordered_lattice = TEAPOT_Ring()
    reorderd_nodes_n =  front_end_nodes + back_end_nodes 
    for n in reorderd_nodes_n:
        reordered_lattice.addNode(n)
    
    return reordered_lattice

def build_injection_region_lattice():
    
    lattice = build_lattice()
    drift_a09 = lattice.getNodeForName("DMCV_A09")
    drift_b01 = lattice.getNodeForName("DMCV_B01")
    #Name:DH_A11B, index:27, type: bend teapot
    #Name:DB23, index:28, type: drift teapot
    #Name:DDMCV_X01 index:51 type: drift teapot
    reordered_lattice = reorganize_lattice(lattice, drift_b01, drift_a09)
    print("method returned lattice starting from exit of A09 to entrance of B01")
    return reordered_lattice

def split_node(node: AccNode, max_part_length: float) -> AccNode:
    if max_part_length is not None and max_part_length > 0.0:
        if node.getLength() > max_part_length:
            node.setnParts(1 + int(node.getLength() / max_part_length))
    return node

def do_bunch(bunch, kineticEnergy):
    #adds bunch initalizes bunch with corresponding parameters and adds 1 particle to the bunch
    m_p = mass_proton
    e_c = charge_electron
    bunch.mass(m_p)
    bunch.charge(e_c)
    bunch.getSyncParticle().kinEnergy(kineticEnergy)
    bunch.addParticle(0,0,0,0,0,0)
    
#----------------------------------------setting up kickers functions -------------------------

def calculate_new_momentum(kinetic_energy):
    p = 0
    p = np.sqrt(((kinetic_energy*10**9)**2 + 2*(kinetic_energy*10**9)*(938*10**6))/(3.0*10**8)**2)
    return p

def calcualte_rigidity(charge, momentum):
    #has to use eV/m/s 
    
    if charge == 1:
    
        rigidity = momentum # p/q -> [eV*s]/[m]*[e]-> [T][m]
        return rigidity
    
    elif charge != 1:
        
        rigidity = momentum/charge
        return rigidity
    
def max_angle_deflection(new_kin_energy, charge, current_old, current_new, selected_plane):
    
    old_max_angles = [12.85*10**-3, 12.84*10**-3, 7.13*10**-3, 7.12*10**-3]
    rigidity_old = 5.657 #T*m from documentation
    alpha_constants = [(rigidity_old/current_old)*x for x in old_max_angles]
    momentum_new = calculate_new_momentum(new_kin_energy)
    rigidity_new = calcualte_rigidity(charge, momentum_new)    
    new_kicks = [(current_new/rigidity_new)*i for i in alpha_constants]   #  
    
    if selected_plane.lower() == "horizontal":
        return [new_kicks[0], new_kicks[2]]
    elif selected_plane.lower() == "vertical":
        return [new_kicks[1], new_kicks[3]]
    elif selected_plane.lower() == "both":
        return new_kicks #[14h, 14v, 23h, 23v]

def get_kickers(lattice, selected_plane):
    
    nodes = []
    if selected_plane.lower() == "both":
       return get_kickers(lattice, "horizontal") + get_kickers(lattice, "vertical")   #h[10,11,12,13] + v{10,11,12,13]}
   
    elif selected_plane.lower() == "horizontal":
        nodes = ["IKICKH_A10", "IKICKH_A11", "IKICKH_A12", "IKICKH_A13" ]   
        
    elif selected_plane.lower() == "vertical":
        nodes = ["IKICKV_A10","IKICKV_A11", "IKICKV_A12", "IKICKV_A13"]
        
    kicker_list = [lattice.getNodeForName(name) for name in nodes]
    return kicker_list

def obtain_scaled_kicker_strengths(lattice, transverse_direction, Hkicks, Vkicks):
    
    direction = transverse_direction
    lattice = lattice
    
    nodes = []
    if direction.lower() == "both":
        return obtain_scaled_kicker_strengths(lattice, "horizontal", Hkicks, Vkicks) + obtain_scaled_kicker_strengths(lattice, "vertical", Hkicks, Vkicks)
    elif direction.lower() == "horizontal":
        nodes = ["IKICKH_A10","IKICKH_A11","IKICKH_A12","IKICKH_A13"]
        strengths = [14.04e-03, -4.28e-03,-4.36727974875e-03,14.092989681e-03]
        angle_kicks = Hkicks        
    elif direction.lower() == "vertical":
        nodes = ["IKICKV_A10","IKICKV_A11","IKICKV_A12","IKICKV_A13"]
        strengths = [8.84e-03, -5.06e-03, -5.32217284098e-03, 9.0098984536e-03]
        angle_kicks = Vkicks
    
    kicks_values = [angle_kicks[0], angle_kicks[1], angle_kicks[1], angle_kicks[0]] #14,23,23,14
    #------ obtain constants or limiting constant values ----------------------------------
    cs = [strengths[i]/kicks_values[i] for i in range(4)] #kickers 1,2,3,4 
    
    # now i have to make sure that all of these are bellow 1 or equal to 1
    #--------------------------------------------------------------------------------------
    bump_array  = np.array(cs) #1,2,3,4
    if np.any(bump_array > 1):
        bump_array /= np.max(bump_array)
    #-----------------------------------------------------
    #obtain new max strengths
    #------------------------------------------------------
    new_strengths = [kicks_values[i] * bump_array[i] for i in range(4)]
    return new_strengths
        
def set_kicker_strengths(lattice, selected_plane, k_strengths):
    
    node = []
    if selected_plane.lower() == "both":
        set_kicker_strengths(lattice, "horizontal", k_strengths[:int(len(k_strengths)/2)])
        set_kicker_strengths(lattice, "vertical", k_strengths[int(len(k_strengths)/2):])
        
    elif selected_plane.lower() == "horizontal":
        node = get_kickers(lattice, selected_plane)
        [node.setParam("kx", strength) for node, strength in zip(node, k_strengths)] #zip pairs 2 entries in lists, (node, k_strength) etc.

    elif selected_plane.lower() == "vertical":
        node = get_kickers(lattice, selected_plane)
        [node.setParam("ky", strength) for node, strength in zip(node, k_strengths)]
#-----------------diagnostics---------------------------------------------------------------------
class BPMSignal:
    """
    This class delivers the average value for coordinate x and y
    """

    def __init__(self):
        self.bunchtwissanalysis = BunchTwissAnalysis()
        self.xAvg = 0.0
        self.yAvg = 0.0
        self.xpAvg = 0.0
        self.ypAvg = 0.0

    def analyzeSignal(self, bunch):
        self.bunchtwissanalysis.analyzeBunch(bunch)

        # if mpi operations are enabled, this section of code will
        # determine the rank of the present node
        rank = 0  # default is primary node
        mpi_init = orbit_mpi.MPI_Initialized()
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        if mpi_init:
            rank = orbit_mpi.MPI_Comm_rank(comm)

        # only the primary node needs to output the calculated information
        if rank == 0:
            self.xAvg = self.bunchtwissanalysis.getAverage(0)
            self.xpAvg = self.bunchtwissanalysis.getAverage(1)
            self.yAvg = self.bunchtwissanalysis.getAverage(2)
            self.ypAvg = self.bunchtwissanalysis.getAverage(3)

    def getSignalX(self):
        return self.xAvg

    def getSignalXP(self):
        return self.xpAvg

    def getSignalY(self):
        return self.yAvg

    def getSignalYP(self):
        return self.ypAvg

class TeapotBPMSignalNode(DriftTEAPOT):
    def __init__(self, name="BPMSignal no name"):
        """
        Constructor. Creates the StatLats TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.bpm = BPMSignal()
        self.setType("BPMSignal")
        self.lattlength = 0.0
        self.setLength(0.0)
        self.position = 0.0

    def track(self, paramsDict):
        """
        The bunchtuneanalysis-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        node = paramsDict["node"]
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.bpm.analyzeSignal(bunch)

    def setPosition(self, pos):
        self.position = pos

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength

    def getSignal(self):
        xAvg = self.bpm.getSignalX()
        yAvg = self.bpm.getSignalY()
        xpAvg = self.bpm.getSignalXP()
        ypAvg = self.bpm.getSignalYP()
        
        return xAvg, xpAvg, yAvg, ypAvg
    
def phase_space_coordinates(lattice, bpm):
    
    for node in lattice.getNodes():
        if node.getName() == bpm:
            return [node.bpm.getSignalX(),node.bpm.getSignalXP(),node.bpm.getSignalY(),node.bpm.getSignalYP()]

def place_diagnostic_node(node, name_of_node, place):
    
    if place.lower() == "entrance":
        bpm_diagnostic = TeapotBPMSignalNode(name_of_node)
        node.addChildNode(bpm_diagnostic, AccNode.ENTRANCE)
    elif place.lower() == "body":
        bpm_diagnostic = TeapotBPMSignalNode(name_of_node)
        node.addChildNode(bpm_diagnostic, AccNode.BODY)
    elif place.lower() == "exit":
        bpm_diagnostic = TeapotBPMSignalNode(name_of_node)
        node.addChildNode(bpm_diagnostic, AccNode.EXIT)

def determine_distance(nodes, target_node):
    pos = 0
    for node in nodes:
        if node.getName() == target_node.getName():
            return pos
        
        pos += node.getLength()
    return None

def get_phase_space_at_diagnostic(lattice, nodes, target_node, diagnostic_node_name, positionalArgument):
    
    diagnostic_node = lattice.getNodeForName(diagnostic_node_name)
    pos = determine_distance(nodes, target_node)
    if positionalArgument.lower() == "entrance":
        children = target_node.getChildNodes(AccNode.ENTRANCE)  
    elif positionalArgument.lower() == "exit": 
        children = target_node.getChildNodes(AccNode.EXIT)
    elif positionalArgument.lower() == "body":
        children = target_node.getChildNodes(AccNode.BODY)   
    for child in children:
        if diagnostic_node_name == child.getName():
            x,xp,y,yp = child.getSignal()
    return x, xp, y, yp, pos

def get_coords(lattice):
    posx, posxp, posy, posyp, pos_z = [], [], [], [], []
    pos = 0.0
    nodes = lattice.getNodes()
    for node in nodes:
        # Check both entrance and exit children AccNode.ENTRANCE, , AccNode.EXIT
        for place in (AccNode.ENTRANCE, AccNode.EXIT):
            for child in node.getChildNodes(place):
                if isinstance(child, TeapotBPMSignalNode):
                    #print(f"currently in {node.getName()} at {place} index {lattice.getNodeIndex(node)}")
                    #print(f"child name: {child.getName()}")
                    x, xp, y, yp = child.getSignal()
                    posx.append(x)
                    posxp.append(xp)
                    posy.append(y)
                    posyp.append(yp)
                    pos_z.append(pos)
        # Always advance lattice position
        pos += node.getLength()
    return [posx, posxp, posy, posyp, pos_z]


def check_lims(k_strengths_old, k_strengths_new):
    #unipolar power source basically means that our kickers can only go in one direction
    success = True

    for i, k in enumerate(k_strengths_new):
        
        this_success=False
        limit = k_strengths_old[i]  
        
        if limit > 0 and 0 < k <= limit:
            this_success = True
        elif limit > 0 and k > limit:
            this_success = False
        elif limit < 0 and limit <= k < 0:
            this_success = True
        elif limit < 0 and k < limit:
            this_success = False
        elif k == 0:
            this_success = False
            
        success = (this_success and success) 
       
    return success

def check_lims_bipolar(k_strengths_old, k_strengths_new):
    
    success = True
    
    for i, k in enumerate(k_strengths_new):
        
        this_success = False
        abs_limit = abs(k_strengths_old[i])
        
        if abs_limit > 0:
            if abs(k) <= abs_limit:
                this_success=True
        elif abs_limit == 0:
            if abs(k) == 0:
                this_success = True
        else:
            this_success = False
            
        success = (this_success and success)
        
    return success

def get_phaseSpaceCoords(lattice, kicker_strengths, target_node, diagnostic_node_name, positionalArgument):
    
    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(1.300)
    bunch.addParticle(0,0,0,0,0,0)
    
    nodes = lattice.getNodes()
    set_kicker_strengths(lattice, "horizontal", kicker_strengths)
    lattice.trackBunch(bunch)
    coords = get_phase_space_at_diagnostic(lattice, nodes, target_node, diagnostic_node_name, positionalArgument)
    return coords

def returnMtrx(matrix):
    
    n = matrix.size()[0]
    m = matrix.size()[1]
    trmatrix = np.zeros([n,m])
    for i in range(n):
        for j in range(m):
            trmatrix[i,j] = matrix.get(i,j)
    
    return trmatrix
    #return matrix
    
def set_all_diagnostics(lattice, selected_plane, target_node):
    #this sets the diagnostics for my problem wont work with any other need a more complex function
    nodes = lattice.getNodes()
    hkickers = get_kickers(lattice, selected_plane)
    for node in nodes:
        if node.getName() == "DMCV_A09":
            place_diagnostic_node(node, f"diagnostic_{node.getName()}", "exit")
        elif node.getName() == "DMCV_B01":
            place_diagnostic_node(node, f"diagnostic_{node.getName()}", "entrance")
        elif node == hkickers[0]:
            place_diagnostic_node(node, f"diagnostic_kicks_{node.getName()}", "entrance")
        elif node in [hkickers[1],hkickers[2]]:
            place_diagnostic_node(node, f"diagnostic_kicks_{node.getName()}", "entrance")
        elif node == hkickers[3]:
            place_diagnostic_node(node, f"diagnostic_kicks_{node.getName()}", "exit")
        elif node == target_node:
            place_diagnostic_node(target_node, f"foil_bpm_{target_node.getName()}", "entrance")
        else:
            place_diagnostic_node(node, f"diagnostic_{node.getName()}", "entrance")
            
#-----------------------------timing fucntionality-----------------------------
@contextmanager
def check_time(label):
    
    start_time = t.time()
    yield
    end_time = t.time()
    elapsed_time = end_time - start_time
    if elapsed_time < 60:
        print(f"[Timer] {label}: {elapsed_time} seconds")
    else:
        print(f"[Timer] {label}: {elapsed_time/60} minutes ")
        
#
#-----------------------------ploting-----------------------------------------
def set_optimizer_plot(lattice, selected_plane, target_node, kicker_strengths, x, xp):
    hkickers = get_kickers(lattice, selected_plane)
    set_all_diagnostics(lattice, selected_plane, target_node)
    new_kicker_strengths = ptf.do_optimize(x, xp, lattice, selected_plane, target_node, kicker_strengths)
    plo.plot_w_new_strengths_wClicker(lattice, hkickers, new_kicker_strengths,[round(x, 6),round(xp, 6)], "October2CoImageClickImages")
    
        
    

# def check_lims(k_strengths_old, k_strengths_new):
#     """
#     Check if all new kicker strengths are within the limits defined
#     by the old kicker strengths.
#     """
#     for i, k in enumerate(k_strengths_new):
#         limit = k_strengths_old[i]

#         if limit > 0:
#             # k must be strictly > 0 and <= limit
#             if k <= 0 or k > limit:
#                 return False

#         elif limit < 0:
#             # k must be strictly < 0 and >= limit
#             if k >= 0 or k < limit:
#                 return False

#         else:
#             # limit == 0 → only valid k is exactly 0
#             if k != 0:
#                 return False

#     return True

            
    
    