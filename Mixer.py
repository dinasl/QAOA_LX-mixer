# import networkx as nx
import numpy as np
# from dataclasses import dataclass, field
from itertools import permutations
# from collections import defaultdict
from typing import Dict, List, Tuple

from Stabilizer import *
from utils import *

# TODO: Implement method for directed graphs (digraph=True). Only for visual representation.
# TODO: Implement method for reduced graphs (reduced=True)
# TODO: Implement blacklist and whitelist methods

class LXMixer:
    def __init__(self, B, nL, digraph=False, reduced=True, sort=False, blacklist=[], whitelist=None):
        self.setB(B, nL, sort)
        
        self.digraph=digraph
        self.reduced=reduced
        self.base_solution_reduced = []
        self.base_solution = []
        self.base_cost=0
        self.base_nedges=0
        
        self.family_of_valid_graphs : Dict[int, List[Tuple[int,...]]] = {}
        self.family_of_valid_graphs_flattened : Dict[int, List[int]] = {}
        self.node_connectors : Dict[int, Dict[int, int]] = {} # Maps nodes to connected nodes to the logical X operators that connect them
        for i in range(self.nB): self.node_connectors[i] = {} # Initializes the node_connectors dictionary for each node

    # Main loop:
        # self.B -> self.compute_family_of_valid_graphs() -> self.family_of_valid_graphs
        # S = Stabilizer(self.B, self.nL, self.restricted)
        # self.family_of_valid_graphs -> S.compute_all_orbits()
            # S : -> self.orbits, self.edges
            #     -> self.num_orbits
        # S.compute_minimal_generating_sets()
            # S : -> self.minimal_generating_sets
        # S.compute_projectors(self.restricted)
            # S : -> self.projectors
        # S.compute_costs()
            # S : -> self.costs
        # best_Gs, best_Xs, best_Zs, best_costs = self.find_best_mixer(S)

    def setB(self, B, nL, sort:bool):
        """
        Sets the feasible set B for the mixer.

        Args:
            B (array-like[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        """
        if isinstance(B, set):
            B = list(B)
        elif isinstance(B, list):
            # B = list(set(B))
            seen = set()
            B = [x for x in B if not (x in seen or seen.add(x))]
        else:
            raise TypeError("B must be a list or a set.")
        self.nB = len(B) # |B|
        if sort:
            B=sorted(B, key=lambda x: int(x, 2))

        if len(B) < 2:
            raise Exception("B must contain at least two elements.")

        self.nL = nL
        for b in B:
            if b >= (1 << self.nL):
                raise Exception(f"Entry {b} exceeds {self.nL} bits.")

        self.B = B
        # Print the binary representation of each element in B as the order may have been reversed
        for b in B:
            print(f"{b:04b}")  
            
    def compute_family_of_valid_graphs(self):
        """
        Computes the family of valid mixers for the feasible set B.
        
        Args:
            B (list[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        
        Returns:
            Dict[int, List[Tuple[int,...]]]: A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator.
        """
        used_X = set()
        
        for i in range(self.nB):
            for j in range(i+1, self.nB):
                X_ij = self.B[i] ^ self.B[j]
                if X_ij not in used_X:
                    self.family_of_valid_graphs[X_ij] = [(i,j)]
                    self.family_of_valid_graphs_flattened[X_ij] = [i,j]
                else: 
                    self.family_of_valid_graphs[X_ij].append((i,j))
                    self.family_of_valid_graphs_flattened[X_ij].extend([i,j])
                used_X.add(X_ij)
        
    def compute_all_orbits(self): # when in the Stabilizer class, args should be self.family_of_valid_graphs
    # DEPTH-FIRST SEARCH ALGORITHM
        
        self.orbits = []
        self.nodes = []
                
        # Maps for each node the |B|-1 logical X operators that connects it to the other nodes
        for X, E in self.family_of_valid_graphs.items():
            for i, j in E:
                self.node_connectors[i][j] = X
                self.node_connectors[j][i] = X

        processed_nodes = set()
        processed_prefixes = set()
        # failed_prefixes = set() # Set to store prefixes of failed paths to avoid using these combinations again

        for seed in range(self.nB):
            if seed in processed_nodes: # Check if the seed node has already been processed
                continue
            
            print(f"\nProcessing seed node {seed}")
            
            seed_Xs = list(self.node_connectors[seed].values()) # All X operators connected to the seed node
            
            stack = []
            initial_available = seed_Xs
            stack.append(([], initial_available, set([seed]))) # Initializes depth-first search with an empty path, available X operators and the seed node
            
            while stack:
                current_path, available_Xs, current_nodes = stack.pop() # Processes each state in last-in-first-out order
                # print("\nCurrent path:")
                # for X in current_path:
                #     print(f"{X:{self.nL}b}")
                # print(f"\nAvailable Xs:")
                # for X in available_Xs:
                #     print(f"{X:0{self.nL}b}")
                # print(f"\nCurrent nodes: {current_nodes}")
                # current_path: sequence of X operators applied so far
                # available_Xs: set of X operators that can still be applied
                # current_nodes: current nodes in orbit
                                
                path_tuple = tuple(sorted(current_path))
                for X in path_tuple: print(f"{X:0{self.nL}b}", end=" ")
                
                if path_tuple in processed_prefixes:
                    # print("if")
                    continue #Continue works?
                processed_prefixes.add(path_tuple) # Add the current path to the processed prefixes
                processed_nodes.update(current_nodes) # Update the processed nodes with the current nodes
                # print(processed_nodes)
                
                if (len(current_nodes) > 2 and len(current_path) > 1): # If the current path has more than one X operator and the current nodes are more than 2, it is a valid orbit
                    if not any(path_tuple == tuple(sorted(orbit)) for orbit in self.orbits): # Check if the orbit is already recorded
                        self.orbits.append(list(current_path))
                        self.nodes.append(list(current_nodes))
                
                if len(current_path) == len(seed_Xs): # If all available X operators have been used, go to a new branch
                    continue
                
                for x, X in enumerate(available_Xs):
                    new_path = current_path + [X]
                    new_available = available_Xs[x+1:] # This is too strict available_Xs - {X}
                    
                    if len(current_path) == 0:
                        new_nodes = self.family_of_valid_graphs_flattened[X] # If no path has been taken yet, the new nodes are the ones connected by the first X operator
                    else: 
                        new_nodes = list({node for u, v in self.family_of_valid_graphs[X] if u in current_nodes and v in current_nodes for node in (u, v)})
                        
                        if not new_nodes:
                            continue
                                            
                    # if len(new_nodes) < 2: # If the new nodes are less than 2, the path is not valid
                    #     continue
                    
                    # processed_prefixes.add(tuple(sorted(new_path)))
                    stack.append((new_path, new_available, new_nodes)) # Add the new path to the stack
                
            
     
# Test family of valid graphs computation
B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
# B = [
#     0b00001,  # weight 1
#     0b00010,
#     0b00100,
#     0b01000,
#     0b10000,
#     0b00011,  # weight 2
#     0b00101,
#     0b00110,
#     0b01001,
#     0b01010,
#     0b01100,
#     0b10001,
#     0b10010,
#     0b10100,
#     0b11000
# ]
# B = [
#     0b10011,
#     0b01100,
#     0b11000,
#     0b00011,
#     0b01001,
#     0b10100,
#     0b00110,
#     0b01110
# ]
lxmixer = LXMixer(B, 4)

lxmixer.compute_family_of_valid_graphs()

print("\nFamily of valid graphs:")
for k, v in lxmixer.family_of_valid_graphs.items():
    print(f"{k:0{lxmixer.nL}b} : {v}")
    
# Test orbit computation
lxmixer.compute_all_orbits()

print("\nnode_connectors:")
for k, v in lxmixer.node_connectors.items():
    print(f"{k}")
    for neighbor, X in v.items():
        print(f"  <-> {neighbor} {X:0{lxmixer.nL}b}")

print("Orbits:")
for orbit in lxmixer.orbits:
    print([f"{X:0{lxmixer.nL}b}" for X in orbit])
print("Orbit nodes:")
print(lxmixer.nodes)