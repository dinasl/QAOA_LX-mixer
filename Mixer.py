# import networkx as nx
import numpy as np
from typing import Dict, List, Tuple
from dataclasses import dataclass

from Stabilizer import *
from utils import *

# TODO: Implement method for directed graphs (digraph=True). Only for visual representation.
# TODO: Implement method for reduced graphs (reduced=True)
# TODO: Implement blacklist and whitelist methods
# TODO: Put the compute_all_orbits method into the Stabilizer class
# TODO: Implement find_best_mixer method

@dataclass
class Orbit:
    """
    Class to store orbits and their properties.
    
    Attributes:
        Xs (List[int]): Logical X operators.
        X_costs List[int]: Logical X operator costs.
        Zs (List[tuple[int, int]]): Projectors (Z operators).
        Z_cost (int): Total cost of the projectors.
    """
    Xs: List[int]  = None
    X_costs: List[int] = None
    Zs: List[Tuple[int, int]] = None
    Z_cost: int = None

class LXMixer:
    """
    Logical X mixer for the QAOA.
    
    The mixer is based on logical X operators connecting nodes in a feasible set B, the span of which is the feasible solution space
    of the QAOA problem. The mixer computes the family of valid graphs, orbit subgraphs within this family and finds the minimal cost combinations of
    subgraphs (with their respective logical X operators and projectors) that connect all nodes in the feasible set B.
    
    Attributes:
        B (list[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        nL (int): Number of qubits of each state. 
        digraph (bool): Whether to use directed graphs (default: False).
        reduced (bool): Whether to use reduced graphs (default: True).
        sort (bool): Whether to sort the feasible set B (default: False).
        blacklist (list[int]): List of logical X operators to exclude from the mixer.
        whitelist (list[int]): List of logical X operators to include in the mixer (default: None).
        
        family_of_valid_graphs (Dict[int, List[Tuple[int,...]]]): A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator.
        S (Stabilizer): Stabilizer object that computes orbits, edges, minimal generating sets, projectors and costs.
        
    Methods:
        setB(): Sets the feasible set B for the mixer.
        compute_family_of_valid_graphs(): Computes the family of valid mixers for the feasible set B
        find_best_mixer(): Finds the best mixer based on the computed orbits, edges, minimal generating sets, projectors and costs.

    """
    def __init__(self, B, nL, digraph=False, reduced=True, sort=False, blacklist=[], whitelist=None):
        self.setB(B, nL, sort)
        
        self.digraph=digraph
        self.reduced=reduced
        self.base_solution_reduced = []
        self.base_solution = []
        self.base_cost=0
        self.base_nedges=0
        
        self.family_of_valid_graphs : Dict[int, List[Tuple[int,...]]] = {}
        # self.family_of_valid_graphs_flattened : Dict[int, List[int]] = {}
        self.node_connectors : Dict[int, Dict[int, int]] = {} # Maps nodes to connected nodes to the logical X operators that connect them
        for i in range(self.nB): self.node_connectors[i] = {} # Initializes the node_connectors dictionary for each node

    # Main loop:
        # self.B -> self.compute_family_of_valid_graphs() -> self.family_of_valid_graphs
        # S = Stabilizer(self.B, self.nL)
        # self.family_of_valid_graphs -> S.compute_all_orbits()
            # S : -> self.orbits, self.edges
            #     -> self.num_orbits
        # S.compute_minimal_generating_sets()
            # S : -> self.minimal_generating_sets
        # S.compute_projectors()
            # S : -> self.projectors
        # S.compute_costs()
            # S : -> self.costs
        # best_Xs, best_Zs, best_costs = self.find_best_mixer(S)

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
        # for b in B:
        #     print(f"{b:04b}")  
            
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
                    # self.family_of_valid_graphs_flattened[X_ij] = [i,j]
                else: 
                    self.family_of_valid_graphs[X_ij].append((i,j))
                    # self.family_of_valid_graphs_flattened[X_ij].extend([i,j])
                used_X.add(X_ij)
        
    def compute_all_orbits(self): # When in the Stabilizer class, args should be self.family_of_valid_graphs
        """
        Computes all orbits in the family of valid graphs using a depth-first search algorithm.
        
        Args:
            family_of_valid_graphs (Dict[int, List[Tuple[int,...]]]): A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator.
            nB (int): Number of nodes in the feasible set B.
            
        """
    # DEPTH-FIRST SEARCH ALGORITHM
        
        self.orbits : Dict[List[int], set[int]] = {}
        # self.orbits = []
        # self.nodes = []
                
        # Maps for each node the |B|-1 logical X operators that connects it to the other nodes
        for X, E in self.family_of_valid_graphs.items():
            for i, j in E:
                self.node_connectors[i][j] = X
                self.node_connectors[j][i] = X

        processed_nodes = set() # Nodes in (>2) orbit 
        processed_prefixes = set()

        for seed in range(self.nB):
            if seed in processed_nodes: # Skip to the next seed that is not already in a (>2) orbit
                continue
            
            seed_Xs = list(self.node_connectors[seed].values()) # All |B|-1 logical X operators connected to the seed node
            
            stack = [] # Stack for depth-first search
            stack.append(([], seed_Xs, set([seed]))) # Initialize
            
            while stack:
                current_path, available_Xs, current_nodes = stack.pop() # Bakctracking step: processes each state in last-in-first-out order
                # current_path: sequence of X operators applied so far
                # available_Xs: set of X operators that can still be applied
                # current_nodes: current nodes in orbit
                                
                path_tuple = tuple(sorted(current_path))
                current_nodes_tuple = tuple(sorted(current_nodes))
                
                if path_tuple in processed_prefixes:
                    continue # If path has already been processed, backtrack
                processed_prefixes.add(path_tuple) # Add the current path to the processed prefixes
                
                if (len(current_nodes) > 2 and len(current_path) > 1): # If the current path has more than one X operator and the current nodes are more than 2, it is a valid orbit
                    # if not any(path_tuple == tuple(sorted(orbit)) for orbit in self.orbits): # Check if the orbit is already recorded
                    #     self.orbits.append(list(current_path))
                    #     self.nodes.append(list(current_nodes))
                    
                    if not any(path_tuple == tuple(sorted(orbit)) for orbit in self.orbits.values()):
                        current_nodes_tuple = tuple(sorted(current_nodes))  # Convert current_nodes to a tuple
                        
                        if current_nodes_tuple in self.orbits.keys():
                            self.orbits[current_nodes_tuple].extend(current_path)  # Extend the list of operators
                            self.orbits[current_nodes_tuple] = sorted(set(self.orbits[current_nodes_tuple]))  # Remove duplicates and sort
                        else:
                            self.orbits[current_nodes_tuple] = sorted(set(current_path))  # Initialize as a sorted list without duplicates
                        
                        processed_nodes.update(current_nodes)  # Update the processed nodes with the current nodes
                
                if len(current_path) == len(seed_Xs): # If the path is |B|-1 long, backtrack
                    continue
                
                for x, X in enumerate(available_Xs): # Iterate over all the next paths in decision tree
                    new_path = current_path + [X]
                    new_available = available_Xs[x+1:]
                    
                    if len(current_path) == 0: # If no path has been taken yet, the new nodes are the ones connected by the first X operator
                        new_nodes = [node for u, v in self.family_of_valid_graphs[X] for node in (u, v)]                    
                    else: 
                        # Sets of nodes that form orbits cannot have edges going out of them
                        new_nodes = list({node for u, v in self.family_of_valid_graphs[X] if u in current_nodes and v in current_nodes for node in (u, v)})
                        
                        if not new_nodes: # If the path doesn't lead anywhere, don't add it to the stack
                            continue
                    
                    stack.append((new_path, new_available, new_nodes)) # Add valid paths to the stack
            
            # If the seed node is not part of any larger (>2) orbit, add all the |B|-1 trivial orbits connecting it
            if seed not in processed_nodes:
                for neighbor, X in self.node_connectors[seed].items():
                    # self.orbits.append([X])
                    # self.nodes.append(sorted([seed, neighbor]))
                    self.orbits[tuple(sorted([seed, neighbor]))] = [X]
        
    def find_best_mixer(self, S: Stabilizer):
        return

# def compute_all_orbits(self, family_of_valid_graphs: Dict[int, List[Tuple[int,...]]], nB): # Function can be put into Stabilizer class
#     """
#     Computes all orbits in the family of valid graphs using a depth-first search algorithm.
    
#     Args:
#         family_of_valid_graphs (Dict[int, List[Tuple[int,...]]]): A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator.
#         nB (int): Number of nodes in the feasible set B.
        
#     """
    
#     node_connectors = {i: {} for i in range(len(family_of_valid_graphs))}
        
#     self.orbits : Dict[List[int], set[int]] = {}
            
#     # Maps for each node the |B|-1 logical X operators that connects it to the other nodes
#     for X, E in family_of_valid_graphs.items():
#         for i, j in E:
#             node_connectors[i][j] = X
#             node_connectors[j][i] = X

#     processed_nodes = set() # Nodes in (>2) orbit 
#     processed_prefixes = set()

#     for seed in range(nB):
#         if seed in processed_nodes: # Skip to the next seed that is not already in a (>2) orbit
#             continue
        
#         seed_Xs = list(node_connectors[seed].values()) # All |B|-1 logical X operators connected to the seed node
        
#         stack = [] # Stack for depth-first search
#         stack.append(([], seed_Xs, set([seed]))) # Initialize
        
#         while stack:
#             current_path, available_Xs, current_nodes = stack.pop() # Bakctracking step: processes each state in last-in-first-out order
#             # current_path: sequence of X operators applied so far
#             # available_Xs: set of X operators that can still be applied
#             # current_nodes: current nodes in orbit
                            
#             path_tuple = tuple(sorted(current_path))
#             current_nodes_tuple = tuple(sorted(current_nodes))
            
#             if path_tuple in processed_prefixes:
#                 continue # If path has already been processed, backtrack
#             processed_prefixes.add(path_tuple) # Add the current path to the processed prefixes
            
#             if (len(current_nodes) > 2 and len(current_path) > 1): # If the current path has more than one X operator and the current nodes are more than 2, it is a valid orbit
#                 # if not any(path_tuple == tuple(sorted(orbit)) for orbit in self.orbits): # Check if the orbit is already recorded
#                 #     self.orbits.append(list(current_path))
#                 #     self.nodes.append(list(current_nodes))
                
#                 if not any(path_tuple == tuple(sorted(orbit)) for orbit in self.orbits.values()):
#                     current_nodes_tuple = tuple(sorted(current_nodes))  # Convert current_nodes to a tuple
                    
#                     if current_nodes_tuple in self.orbits.keys():
#                         self.orbits[current_nodes_tuple].extend(current_path)  # Extend the list of operators
#                         self.orbits[current_nodes_tuple] = sorted(set(self.orbits[current_nodes_tuple]))  # Remove duplicates and sort
#                     else:
#                         self.orbits[current_nodes_tuple] = sorted(set(current_path))  # Initialize as a sorted list without duplicates
                    
#                     processed_nodes.update(current_nodes)  # Update the processed nodes with the current nodes
            
#             if len(current_path) == len(seed_Xs): # If the path is |B|-1 long, backtrack
#                 continue
            
#             for x, X in enumerate(available_Xs): # Iterate over all the next paths in decision tree
#                 new_path = current_path + [X]
#                 new_available = available_Xs[x+1:]
                
#                 if len(current_path) == 0: # If no path has been taken yet, the new nodes are the ones connected by the first X operator
#                     new_nodes = [node for u, v in family_of_valid_graphs[X] for node in (u, v)]                    
#                 else: 
#                     # Sets of nodes that form orbits cannot have edges going out of them
#                     new_nodes = list({node for u, v in family_of_valid_graphs[X] if u in current_nodes and v in current_nodes for node in (u, v)})
                    
#                     if not new_nodes: # If the path doesn't lead anywhere, don't add it to the stack
#                         continue
                
#                 stack.append((new_path, new_available, new_nodes)) # Add valid paths to the stack
        
#         # If the seed node is not part of any larger (>2) orbit, add all the |B|-1 trivial orbits connecting it
#         if seed not in processed_nodes:
#             for neighbor, X in node_connectors[seed].items():
#                 # self.orbits.append([X])
#                 # self.nodes.append(sorted([seed, neighbor]))
#                 self.orbits[tuple(sorted([seed, neighbor]))] = [X]

# Standalone code (e.g., debugging or testing)
if __name__ == '__main__':
    # B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
    # B = [
    #     0b00001,
    #     0b00010,
    #     0b00100,
    #     0b01000,
    #     0b10000,
    #     0b00011,
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
    B = [
        0b10011,
        0b01100,
        0b11000,
        0b00011,
        0b01001,
        0b10100,
        0b00110,
        0b01110
    ]
    lxmixer = LXMixer(B, 5)

    lxmixer.compute_family_of_valid_graphs()

    print("\nFamily of valid graphs:")
    for k, v in lxmixer.family_of_valid_graphs.items():
        print(f"{k:0{lxmixer.nL}b} : {v}")

    lxmixer.compute_all_orbits()

    print("\nnode_connectors:")
    for k, v in lxmixer.node_connectors.items():
        print(f"{k}")
        for neighbor, X in v.items():
            print(f"  <-> {neighbor} {X:0{lxmixer.nL}b}")

    print("Orbits:")
    # for orbit in lxmixer.orbits:
    #     print([f"{X:0{lxmixer.nL}b}" for X in orbit])
    # print("Orbit nodes:")
    # print(lxmixer.nodes)
    # print(lxmixer.orbits)

    for nodes, Xs in lxmixer.orbits.items():
        print(f"{nodes} : [{', '.join(f'{X:0{lxmixer.nL}b}' for X in Xs)}]")