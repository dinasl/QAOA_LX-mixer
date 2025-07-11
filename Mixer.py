# import networkx as nx
# import numpy as np
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass, field
from itertools import combinations
import math

from Stabilizer import *
from utils import ncnot, is_connected

# TODO: Implement method for directed graphs (digraph=True). Only for visual representation.
# TODO: Implement method for reduced graphs (reduced=True)
# TODO: Implement blacklist and whitelist methods

@dataclass
class Orbit:
    """
    Class to store orbits and their properties.
    
    Attributes:
        Xs (Set[int]): Logical X operators.
        X_costs List[int]: Logical X operator costs.
        Zs (List[tuple[int, int]]): Projectors (Z operators).
        Z_cost (int): Total cost of the projectors.
    """
    Xs: Set[int] = field(default_factory=list)
    X_costs: List[int] = field(default_factory=list)
    Zs: List[Tuple[int, int]] = field(default_factory=list)
    Z_cost: int = float('inf')
    cost:int = float('inf')

class LXMixer:
    """
    Logical X mixer for the QAOA.
    
    The mixer is based on logical X operators connecting nodes in a feasible set B, the span of which is the feasible solution space
    of the QAOA problem. The mixer computes the family of valid graphs, orbit subgraphs within this family and finds the minimal cost combinations of
    subgraphs (with their respective logical X operators and projectors) that connect all nodes in the feasible set B.
    
    Attributes:
        B (list[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        nB (int): Number of elements in the feasible set B.
        nL (int): Number of qubits. 
        digraph (bool): Whether to use directed graphs (default: False).
        reduced (bool): Whether to use reduced graphs (default: True).
        sort (bool): Whether to sort the feasible set B (default: False).
        blacklist (list[int]): List of logical X operators to exclude from the mixer.
        whitelist (list[int]): List of logical X operators to include in the mixer (default: None).
        
        family_of_valid_graphs (Dict[int, List[Tuple[int,...]]]): A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator.
        orbits (Dict[Tuple[int,...], Orbit]): A dictionary mapping node tuples to Orbit objects, containing the logical X operators, their respective costs, the Z operators representing the orbit and their total cost.
        S (Stabilizer): Stabilizer object that computes the orbit's respective projectors.
        
    Methods:
        setB(): Sets the feasible set B for the mixer.
        compute_family_of_valid_graphs(): Computes the family of valid mixers for the feasible set B
        compute_all_orbits(): Computes all orbits in the family of valid graphs using a depth-first search algorithm.
        compute_costs(): Computes and updates the costs in the Orbit objects in orbits.
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
        # LX = LXMixer(B, nL)
        # LX.compute_family_of_valid_graphs() -> LX.family_of_valid_graphs
        # LX.compute_all_orbits() -> LX.orbits
        # S = Stabilizer(LX.B, LX.nL, LX.orbits)
        # S.compute_minimal_generating_sets() -> updates S.orbits and LX.orbits
        # S.compute_projectors() -> -"-
        # LX.compute_costs() -> -"-
        # best_Xs, best_Zs, best_cost = LX.find_best_mixer()

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
        self.orbits : Dict[Tuple[int,...], Orbit] = {}
        
        # self.orbits : Dict[List[int], set[int]] = {}
        
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
                    
                    if not any(path_tuple == tuple(sorted([X])) for orbit in self.orbits.values() for X in orbit.Xs):
                        # self.orbits.append(list(current_path))
                        # self.nodes.append(list(current_nodes))
                        
                        # Create a new Orbit object if it doesn't exist
                        if current_nodes_tuple not in self.orbits.keys():
                            self.orbits[current_nodes_tuple] = Orbit(Xs=list(current_path))
                        else:
                            self.orbits[current_nodes_tuple].Xs.extend(current_path)
                            self.orbits[current_nodes_tuple].Xs = list(set(self.orbits[current_nodes_tuple].Xs))
                            # self.orbits[current_nodes_tuple].Xs = sorted(set(self.orbits[current_nodes_tuple].Xs)) 
                    
                    # if not any(path_tuple == tuple(sorted(orbit)) for orbit in self.orbits.values()):
                    #     current_nodes_tuple = tuple(sorted(current_nodes))  # Convert current_nodes to a tuple
                        
                    #     if current_nodes_tuple in self.orbits.keys():
                    #         self.orbits[current_nodes_tuple].extend(current_path)  # Extend the list of operators
                    #         self.orbits[current_nodes_tuple] = sorted(set(self.orbits[current_nodes_tuple]))  # Remove duplicates and sort
                    #     else:
                    #         self.orbits[current_nodes_tuple] = sorted(set(current_path))  # Initialize as a sorted list without duplicates
                            
                    # if not any(path_tuple == tuple(sorted(orbit)) for orbit in self.orbits): # Check if the orbit is already recorded
                    #     self.orbits.append(list(current_path))
                    #     self.nodes.append(list(current_nodes))
                        
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
                    self.orbits[tuple(sorted([seed, neighbor]))] = Orbit(Xs=[X])
                    
                    # self.orbits[tuple(sorted([seed, neighbor]))] = [X]
                    
                    # self.orbits.append([X])
                    # self.nodes.append(sorted([seed, neighbor]))
                    
    def compute_costs(self):
        """
        Computes and updates the costs in the Orbit objects in orbits.
        """
        for nodes, orbit in self.orbits.items():
            
            orbit.Z_cost = sum([ncnot(Z[1]) for Z in orbit.Zs]) if orbit.Zs else 0
            
            if math.log2(len(nodes)) < len(orbit.Xs):
                orbit.Xs, orbit.X_costs = zip(*sorted(zip(orbit.Xs, [ncnot(X) for X in orbit.Xs]), key=lambda x: x[1])[:int(math.log2(len(nodes)))])
            else:
                orbit.X_costs = [ncnot(X) for X in orbit.Xs]
            
            orbit.cost = sum(orbit.X_costs) + orbit.Z_cost
        
    def find_best_mixer(self):
        
        best_cost = float('inf')
        best_combinations = []
        
        if len(self.orbits.keys()) == 1:
            best_Xs = [list(self.orbits.values())[0].Xs]
            best_Zs = [list(self.orbits.values())[0].Zs]
            best_cost = sum(list(self.orbits.values())[0].X_costs) # There is no projector needed
            return best_Xs, best_Zs, best_cost
        
        N = range(2, len(self.orbits.keys()))
        for n in N:
            for combination in combinations(self.orbits.keys(), n):
                # print(f"\nChecking combination: {combination}")
                if len(set([node for nodes in combination for node in nodes])) != self.nB:
                    # print(f"Combination does not cover all nodes, skipping.")
                    continue
                if not is_connected(combination):
                    # print(f"Combination is not connected, skipping.")
                    continue
                # print(f"Combination is connected, computing cost...")
                cost = 0
                for orbit_nodes in combination:
                    # cost += sum(self.orbits[orbit_nodes].X_costs[:int(math.log2(len(orbit_nodes)))]) + self.orbits[orbit_nodes].Z_cost
                    cost += self.orbits[orbit_nodes].cost
                    if cost > best_cost:
                        break
                # print(f"Cost: {cost}")    
                if cost < best_cost:
                    # print(f"New best combination found: {combination} with cost {cost}")
                    best_cost = cost
                    best_combinations = [combination]
                elif cost == best_cost:
                    # print(f"Combination {combination} has the same cost as the best one, adding to the list.")
                    best_combinations.append(combination)
        
        best_Xs = [[self.orbits[orbit_nodes].Xs for orbit_nodes in combination] for combination in best_combinations]
        best_Zs = [[self.orbits[orbit_nodes].Zs for orbit_nodes in combination] for combination in best_combinations]
        
        return best_Xs, best_Zs, best_cost    
                    
# def is_connected(orbits):
#     # Short circuit evaluation
#     for orbit in orbits:
#         if not any(set(orbit.intersection(set(other_orbit))) for other_orbit in orbits if other_orbit != orbit):
#             return False
#     return True

# Standalone code

if __name__ == '__main__':
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
    import time
    
    B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
    print(f"\nB = {[f'{b:0{len(bin(max(B)))-2}b}' for b in B]}")
    
    lxmixer = LXMixer(B, 4)
    # lxmixer = LXMixer(B, 5)

    print("\nComputing family of valid graphs...")
    start_time = time.time()
    lxmixer.compute_family_of_valid_graphs()
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")

    print("\nFamily of valid graphs:")
    for k, v in lxmixer.family_of_valid_graphs.items():
        print(f"{k:0{lxmixer.nL}b} : {v}")

    print("\nComputing all orbits...")
    start_time = time.time()
    lxmixer.compute_all_orbits()
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")
    
    # print("\nnode_connectors =")
    # for k, v in lxmixer.node_connectors.items():
    #     print(f"{k}")
    #     for neighbor, X in v.items():
    #         print(f"  <-> {neighbor} {X:0{lxmixer.nL}b}")

    print("\nOrbits (without projectors and costs):")
    for nodes, orbit in lxmixer.orbits.items():
        print(f"{nodes} : [{', '.join(f'{X:0{lxmixer.nL}b}' for X in orbit.Xs)}]")
    
    S = Stabilizer(lxmixer.B, lxmixer.nL, lxmixer.orbits)
    
    print("\nComputing minimal generating sets...")
    start_time = time.time()
    S.compute_minimal_generating_sets()
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")
    
    print("\nComputing projectors...")
    start_time = time.time()
    S.compute_projector_stabilizers()
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")
    
    print("\nComputing costs...")
    start_time = time.time()
    lxmixer.compute_costs()
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")
    
    print("\nOrbits with projectors and costs:")
    for nodes, orbit in lxmixer.orbits.items():
        print(f"{nodes} : Xs = [{', '.join(f'{X:0{lxmixer.nL}b}' for X in orbit.Xs)}], Zs = [{', '.join(f'{"+" if Z[0] == 1 else "-"}{Z[1]:0{lxmixer.nL}b}' for Z in orbit.Zs if len(Z) == 2)}]")
        print(f"cost = {orbit.cost}")
    
    print("\nFinding best mixer...")
    start_time = time.time()
    best_Xs, best_Zs, best_cost = lxmixer.find_best_mixer()
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")
    
    print(f"\nFound {len(best_Xs)} best combinations of orbits with cost {best_cost}.")
    print("\nBest mixer:")
    print(f"[{', '.join(f'[{", ".join(f"[{", ".join(f"{x:0{lxmixer.nL}b}" for x in sub_Xs)}]" for sub_Xs in Xs)}]' for Xs in best_Xs)}]")
    print("\nBest projectors:")
    print(f"[{', '.join(f'[{", ".join(f'{"+" if z[0] > 0 else "-"}{z[1]:0{lxmixer.nL}b}' for z in Z)}]' for sub_Zs in best_Zs for Z in sub_Zs)}]")    