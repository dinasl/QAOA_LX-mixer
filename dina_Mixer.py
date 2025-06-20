import itertools
import networkx as nx
import numpy as np
import json
# from scipy.special import comb
import sys

from PauliString import *
from PauliOperations import *
# from Stabilizers import *
# from GroupGraph import *

import math
#import openquantumcomputing._rust

from tqdm import tqdm

from sympy import *
from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product
#from sympy.physics.quantum import TensorProduct

#from networkx.algorithms.approximation import clique

class Mixer:
    def __init__(self, B, digraph=False, reduced=True, sort=False, blacklist=[], whitelist=None):
        self.setB(B, sort)
        self.digraph=digraph
        self.reduced=reduced
        self.base_solution_reduced = []
        self.base_solution = []
        self.base_cost=0
        self.base_nedges=0
        self.compute_commuting_pairs()

        # self.compute_family_of_graphs(blacklist=blacklist)

        # Graph stuff

    def setB(self, B, sort):
        if isinstance(B, set):
            B = list(B)
        elif isinstance(B, list):
            B = list(set(B))### to make it unique
        else:
            raise TypeError("B must be a list or a set.")
        self.nB = len(B)  ### number of B's
        if sort:
            B=sorted(B, key=lambda x: int(x, 2))

        if len(B) < 2:
            raise Exception("B must contain at least two elements.")

        self.nL = len(B[0])### number of literals
        for b in B:
            if len(b) != self.nL:
                raise Exception("All entries of B must have the same length.")

        self.B = []
        for b in B:
            self.B.append(BitString(1,b))
            
    def compute_commuting_pairs(self):
        self.commuting_pairs = {}
        for i in range(self.nB):
            for j in range(i + 1, self.nB):
                Xij = Xoperator(self.B[i], self.B[j]).P
                self.commuting_pairs[Xij] = self.commuting_pairs.get(Xij, [])
                self.commuting_pairs[Xij].append([i, j])