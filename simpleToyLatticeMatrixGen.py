#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:36:17 2025

@author: l5g
"""
import sympy 
from sympy import Matrix, pprint, init_printing
from sympy import MatrixSymbol, simplify
import matplotlib.pyplot as plt
import numpy as np

from orbit.core.bunch import Bunch
from orbit.teapot import DriftTEAPOT, QuadTEAPOT, KickTEAPOT, BaseTEAPOT
from orbit.teapot.teapot_matrix_lattice import TEAPOT_MATRIX_Lattice, TEAPOT_Lattice
from orbit.matrix_lattice import BaseMATRIX, MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccNodeBunchTracker, AccActionsContainer


def output_matrix(matrix):
    n = matrix.size()[0]
    m = matrix.size()[1]
    trM = np.zeros([n,m])
    for i in range(n):
        for j in range(m):
            trM[i,j] = matrix.get(i,j)
    return trM


bunch = Bunch()
bunch.mass(.938)
bunch.charge(1)
bunch.getSyncParticle().kinEnergy(1.3)
bunch.addParticle(0,0,0,0,0,0)

lattice = TEAPOT_Lattice("toy lattice")
toy_drift = DriftTEAPOT("test drift")
toy_drift.setLength(4)
toyQuad = QuadTEAPOT("quad toy")
toyQuad.setParam("kq", +1.4)
toyQuad.setLength(1)

lattice.addNode(toy_drift)
lattice.addNode(toyQuad)
lattice.initialize()

print(lattice)
print(toy_drift.getLength())

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)

matrix_nodes = matrix_lattice.getNodes()
for n in matrix_nodes:
    print(f"type:{n.getType()}")
    if isinstance(n, BaseMATRIX) == True:
        trMatrix = output_matrix(n.getMatrix())
        symMatrx= Matrix(trMatrix)
        print(f"\nMatrix at node {n.getName()}")
        pprint(symMatrx)
        print()



