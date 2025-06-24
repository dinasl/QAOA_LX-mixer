import networkx as nx
import numpy as np
from dataclasses import dataclass, field
from itertools import combinations, product

from Stabilizer import *
from utils import *

# import json
# from scipy.special import comb
# import sys

# import math
# import openquantumcomputing._rust

# from sympy import *
# from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product
# from sympy.physics.quantum import TensorProduct

# from networkx.algorithms.approximation import clique

@dataclass
class XGraph:
    """
    Attributes:
        X (int): Mixer mask (binary int representation).
        G (nx.Graph): The graph representing the mixer.
        cost (int): The cost associated with the mixer.
    """
    X : int
    G : nx.Graph
    cost : int = None
    subgraphs = field(default_factory=list)
    projectors = field(default_factory=list) # 2D array
    subgraph_costs = field(default_factory=list)
    
@dataclass
class CombinedXGraph:
    """
    Attributes:
        Xs (list[int]): List of mixer masks (binary int representations).
        G (nx.Graph): The combined graph representing the mixer.
        cost (int): The cost associated with the combined mixer.
    """
    Xs : list[int]
    G : nx.Graph
    projectors = field(default_factory=list) # 2D array
    cost : int = None

# TODO: Implement method for directed graphs (digraph=True)
# TODO: Implement method for reduced graphs (reduced=True)
# TODO: Implement blacklist and whitelist methods

class LXMixer:
    def __init__(self, B, digraph=False, reduced=True, sort=False, blacklist=[], whitelist=None):
        self.setB(B, sort)
        
        self.digraph=digraph
        self.reduced=reduced
        self.base_solution_reduced = []
        self.base_solution = []
        self.base_cost=0
        self.base_nedges=0
        
    # Main loop:
        # family_of_valid_graphs = self.compute_family_of_valid_graphs()
        # self.compute_all_subgraphs(family_of_valid_graphs)
        # connected_combinations = self.check_all_combinations(family_of_valid_graphs, combine_subgraphs=True)
        # self.optimal_mixer = min(connected_combinations, key=lambda x: x.cost)

    def setB(self, B, sort:bool):
        if isinstance(B, set):
            B = list(B)
        elif isinstance(B, list):
            B = list(set(B)) # To make it unique
        else:
            raise TypeError("B must be a list or a set.")
        self.nB = len(B) # |B|
        if sort:
            B=sorted(B, key=lambda x: int(x, 2))

        if len(B) < 2:
            raise Exception("B must contain at least two elements.")

        self.nL = len(str(B[0])) # Number of literals
        for b in B:
            if len(str(b)) != self.nL:
                raise Exception("All entries of B must have the same length.")

        self.B = B
            
    def compute_family_of_valid_graphs(self):
        """
        Computes the family of valid mixers for the feasible set B.
        
        Args:
            B (array-like[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        
        Returns:
            list[MixerGraph]: A list of MixerGraph structs representing the family of valid mixers.
        """
        dim = len(self.B)
        used_X = set()
        graphs = []
        
        for i in range(dim):
            for j in range(i+1, dim):
                edges = set()
                X_ij = self.B[i] ^ self.B[j]  # Logical X operator connecting nodes i and j
                if X_ij in used_X:
                    continue
                unconnected_mask = abs(1 - np.isin(self.B, X_ij^self.B)) # Nodes in B that are not connected by X_ij
                b1 = np.delete(self.B, unconnected_mask)  # Nodes in B that are connected by X_ij
                b2 = b1^X_ij  # Their corresponding connected nodes
                edges = list(zip(b1, b2))
                # edges.add(tuple(sorted(u,v) for u, v in zip(b1, b2)))
                graphs.append(XGraph(X=X_ij, G=nx.Graph(edges), cost=ncnot(X_ij)))
                
                used_X.add(X_ij)
        
        return graphs

    def compute_all_subgraphs(self, xgraph_list:list[XGraph]):
        """
        Computes all subgraphs of a given list of XGraph.
        
        Args:
            graph_list (list[XGraphs]]): List of XGraph structs to divide into subgraphs.
        """
        for xG in xgraph_list:
            edges = list(xG.G.edges)
            for r in range(1, len(edges) + 1):
                if is_power_of_two(r):
                    for subset in combinations(edges, r):
                        subG = nx.Graph()
                        subG.add_edges_from(subset)
                        if check_if_orbit(subG.nodes)[0]:
                            xG.subgraphs.append(subG)
                        
                            # TODO: Restricted = True/False option
                            xG.projectors.append(np.array(compute_restricted_projector_stabilizer(subG.nodes(), self.nL)))
                            xG.subgraph_costs.append(np.sum(ncnot(xG.projectors[-1])))
                        
        # pass-by-reference  
    
    def check_all_combinations(self, graph_list:list[XGraph], N=None, combine_subgraphs=True):
        """
        Checks all combinations of (sub)graphs in the given list of XGraph.
        If n is specified, only combinations of size n are checked.

        Args:
            graph_list (list[XGraph]): List of XGraph structs to check combinations from.
            N (int, optional): Size of combinations. Defaults to None.

        Returns:
            list[XGraph]: List of XGraph structs that are connected.
        """
        connected_combinations = []
        
        if N is None: N = range(1,len(graph_list)+1)
        else: N = range(1,N+1)
        
        options = []
        for xG in graph_list:
            if combine_subgraphs:
                options.append([(xG.X, subgraph, projectors_i, subgraph_cost) for subgraph, projectors_i, subgraph_cost in zip(xG.subgraphs, xG.projectors, xG.subgraph_costs)] + [(None, nx.Graph(), None, None)])
            else: options.append([(xG.X, xG.G, None, xG.cost), (None, nx.Graph(), None, None)])
        
        for n in N:
            for subset in combinations(options,n): # Unordered combinations of size n
                for combination in product(*subset):
                    combined_G = nx.Graph()
                    Xs = []
                    projectors = []
                    combined_cost = 0
                    for X_i, subG, projectors_i, cost in combination:
                        combined_G = nx.compose(combined_G, subG)
                        # TODO: Check connected
                        if X_i is not None: Xs.append(X_i)
                        projectors.append(projectors_i)
                        combined_cost += cost
                                                
                    if combined_G.number_of_nodes == len(self.B) and combined_G.number_of_edges == (len(self.B)-1): # Easy check
                        if nx.is_connected(combined_G): # Actual check
                            connected_combinations.append(
                                CombinedXGraph(
                                    Xs=Xs, 
                                    G=combined_G, 
                                    projectors=projectors,
                                    cost = combined_cost
                                    )
                                )
        return connected_combinations