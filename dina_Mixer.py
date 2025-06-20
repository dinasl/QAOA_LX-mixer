import itertools
import networkx as nx
import numpy as np
import json
# from scipy.special import comb
import sys
from dataclasses import dataclass

import math
#import openquantumcomputing._rust

from sympy import *
from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product
#from sympy.physics.quantum import TensorProduct

#from networkx.algorithms.approximation import clique

#
@dataclass
class MixerGraph:
    
    X : int
    G : nx.Graph
    cost : int
    
    """
    Attributes:
        X (int): Mixer mask (integer representation).
        G (nx.Graph): The graph representing the mixer.
        cost (int): The cost associated with the mixer.
    """

class LXMixer:
    def __init__(self, B, digraph=False, reduced=True, sort=False, blacklist=[], whitelist=None):
        self.setB(B, sort)
        self.digraph=digraph
        self.reduced=reduced
        self.base_solution_reduced = []
        self.base_solution = []
        self.base_cost=0
        self.base_nedges=0

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
            
    def compute_family_of_valid_mixergraphs(self):
        """
        Computes the family of valid mixers for the feasible set B.
        
        List of MixerGraph structs. 
        
        Avoid making duplicates (make struct in one go)? If not, set(). 
        """
        pass
    
#def main loop
    # compute family of valid mixer graphs
    #if choose_n: b_func
        #n = 1
        #connected = False
        #while not connected:
            #all combinations of n graphs
            #check if any connected:
                #append to list, connected = True
            #n+=1
    #compute all subgraphs MixerGraphs (with projectors?) compute_all_subgraphs()
    #check all connected combinations with cost
    #take lowest --> output MixerGraph with lowest cost
    
    #projector

def compute_all_subgraphs(mixergraph_list):
    
    """
    Args:
        graph_list (list[MixerGraphs]]): List of MixerGraphs structs.
    """
    pass

def cost():
    pass

def check_if_orbit():
    pass

def algorithm_1():
    pass

def linalg_approach():
    pass

def stabilizer_approach():
    pass