# import networkx as nx
# import numpy as np
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass, field
from itertools import combinations
import math
from tqdm import tqdm
import sys

from Stabilizer import *
from utils import ncnot, is_connected, split_into_suborbits
from plot_mixers import plot_mixer_graph

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
        orbits (Dict[Tuple[int,...], Orbit]): A dictionary mapping node tuples to Orbit objects, containing the logical X operators, the Z operators representing the orbit and their total cost.
        S (Stabilizer): Stabilizer object that computes the orbit's respective projectors.
        
        best_Xs (List[List[List[int]]]): List(s) of lists of logical X operators that form the best mixers.
        best_Zs (List[List[List[Tuple[int, int]]]]): List(s) of lists of projectors (Z operators) that form the best mixers.
        best_combinations (List[List[Tuple[int,...]]]): List(s) of tuples of node indices that form the best mixers.
        best_cost (int): The cost of the best mixer(s) found.
        
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
        self.orbits : Dict[Tuple[int,...], Orbit] = {}
        
        self.best_Xs = []
        self.best_Zs = []
        self.best_combinations = []
        self.best_cost: int = float('inf')

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
        Results stored in the `orbits` attribute as a dictionary mapping node tuples to Orbit objects.
        """    
                
        # Maps for each node the |B|-1 logical X operators that connects it to the other nodes
        for X, E in self.family_of_valid_graphs.items():
            for i, j in E:
                self.node_connectors[i][j] = X
                self.node_connectors[j][i] = X

        processed_nodes = [] # Nodes in (>2) orbit 
        # processed_nodes = set()

        for seed in range(self.nB):
            temp_processed_nodes = []
            
            seed_Xs = list(self.node_connectors[seed].values()) # All |B|-1 logical X operators connected to the seed node
            
            stack = [] # Stack for depth-first search
            stack.append(([], seed_Xs, set([seed]))) # Initialize
            while stack:
                current_path, available_Xs, current_nodes = stack.pop() # Bakctracking step: processes each state in last-in-first-out order
                current_nodes_tuple = tuple(sorted(current_nodes))
                
                if (len(current_nodes) >= 2 and len(current_path) > 1): # If the current path has more than one X operator and the current nodes are more than 2, it is a valid orbit                     
                    
                    # Create a new Orbit object if it doesn't exist
                    if current_nodes_tuple not in self.orbits.keys():
                        self.orbits[current_nodes_tuple] = Orbit(Xs=list(current_path))
                        temp_processed_nodes.append(current_nodes_tuple)
                    else:
                        self.orbits[current_nodes_tuple].Xs.extend(current_path)
                        self.orbits[current_nodes_tuple].Xs = list(set(self.orbits[current_nodes_tuple].Xs))
                        # Update the processed nodes with the current nodes
                
                if len(current_path) == len(seed_Xs): # If the path is |B|-1 long, backtrack
                    continue
                                
                for x, X in enumerate(available_Xs): # Iterate over all the next paths in decision tree
                    new_path = current_path + [X]
                    
                    new_available = available_Xs[:x] + available_Xs[x+1:] # available_Xs - {X}
                    
                    if len(current_path) == 0: # If no path has been taken yet, the new nodes are the ones connected by the first X operator
                        new_nodes = [node for u, v in self.family_of_valid_graphs[X] for node in (u, v)]                    
                    else: 
                        # Sets of nodes that form orbits cannot have edges going out of them
                        new_nodes = [node for u, v in self.family_of_valid_graphs[X] if u in current_nodes and v in current_nodes for node in (u, v)]
                        
                        if len(new_nodes) == 0 or seed not in new_nodes: # If the path doesn't lead anywhere, don't add it to the stack # THIS IS NEWWWW
                            continue
                        
                        if tuple(sorted(new_nodes)) in processed_nodes: # If the new nodes are already processed, skip it
                            continue
                        
                    stack.append((new_path, new_available, new_nodes)) # Add valid paths to the stack
            
            # Split unconnected graphs into suborbits
            for nodes in temp_processed_nodes:
                Xs = self.orbits[nodes].Xs
                if len(Xs) < len(nodes)-1: # If the orbit is not complete
                    new_orbits = split_into_suborbits(self.family_of_valid_graphs, Xs, nodes)
                    self.orbits.pop(nodes) # Remove the orbit from the dictionary
                    for new_orbit in new_orbits:
                        self.orbits[tuple(sorted(new_orbit))] = Orbit(Xs=Xs)

            processed_nodes = list(self.orbits.keys())
            
            if not any(seed in nodes for nodes in processed_nodes):
                for neighbour, X in self.node_connectors[seed].items():
                    self.orbits[tuple(sorted([seed, neighbour]))] = Orbit(Xs=[X])
                    
        print(len(self.orbits.keys()), "orbits. ", len(set(self.orbits.keys())), "unique orbits found.")
        
        # Split unconnected graphs into suborbits
        # orbit_nodes_to_check = list(self.orbits.keys())
        # for nodes in orbit_nodes_to_check:
        #     Xs = self.orbits[nodes].Xs
        #     if len(Xs) < len(nodes)-1:
        #         new_orbits = split_into_suborbits(self.family_of_valid_graphs, Xs, nodes)
        #         self.orbits.pop(nodes)
        #         for new_orbit in new_orbits:
        #             self.orbits[tuple(sorted(new_orbit))] = Orbit(Xs=Xs)
        
        # print(len(self.orbits), same_nodes)
                    
    def compute_costs(self):
        """
        Computes and updates the costs in the Orbit objects in `orbits`.
        """
        for nodes, orbit in self.orbits.items():
            
            orbit.Z_cost = sum([ncnot(Z[1]) for Z in orbit.Zs]) if orbit.Zs else 0
            
        #     # if math.log2(len(nodes)) < len(orbit.Xs):
            orbit.Xs, orbit.X_costs = zip(*sorted(zip(orbit.Xs, [ncnot(X) for X in orbit.Xs]), key=lambda x: x[1])[:int(math.log2(len(nodes)))])
        #     # else:
        #     #     orbit.X_costs = [ncnot(X) for X in orbit.Xs]
            
            orbit.cost = sum(orbit.X_costs) + orbit.Z_cost
            
        # for nodes in self.orbits.keys():
        #     self.orbits[nodes].Z_cost = sum([ncnot(Z[1]) for Z in self.orbits[nodes].Zs]) if self.orbits[nodes].Zs else 0

        #     self.orbits[nodes].Xs, self.orbits[nodes].X_costs = zip(*sorted(zip(self.orbits[nodes].Xs, [ncnot(X) for X in self.orbits[nodes].Xs]), key=lambda x: x[1])[:int(math.log2(len(nodes)))])
            
        #     self.orbits[nodes].cost = sum(self.orbits[nodes].X_costs) + self.orbits[nodes].Z_cost
        
    def find_best_mixer(self):
        """ 
        Finds the best mixer based on the computed orbits, edges, minimal generating sets, projectors and costs.
        """
        
        # best_cost = float('inf')
        # best_combinations = []
        
        if len(self.orbits.keys()) == 1:
            self.best_Xs = [[list(self.orbits.values())[0].Xs]]
            self.best_Zs = [[[(1,0)]]] #TODO
            self.best_cost = sum(list(self.orbits.values())[0].X_costs) # There is no projector needed
            self.best_combinations = [[tuple(self.B)]]
            # return best_Xs, best_Zs, best_cost
            return
        
        N = range(2, len(self.orbits.keys())+1)
        for n in N:
            for combination in combinations(self.orbits.keys(), n):
                # time.sleep(0.05)
                if len(set([node for nodes in combination for node in nodes])) != self.nB:
                    print(f"Combination {combination} does not cover all nodes in B, skipping.")
                    continue
                if not is_connected(combination):
                    continue
                cost = 0
                for orbit_nodes in combination:
                    cost += self.orbits[orbit_nodes].cost
                    if cost > self.best_cost:
                        break
                # print(f"Cost: {cost}")    
                if cost < self.best_cost:
                    # print(f"New best combination found: {combination} with cost {cost}")
                    self.best_cost = cost
                    self.best_combinations = [combination]
                elif cost == self.best_cost:
                    # print(f"Combination {combination} has the same cost as the best one, adding to the list.")
                    self.best_combinations.append(combination)
        
        self.best_Xs = [[self.orbits[orbit_nodes].Xs for orbit_nodes in combination] for combination in self.best_combinations]
        self.best_Zs = [[self.orbits[orbit_nodes].Zs for orbit_nodes in combination] for combination in self.best_combinations]
        
        return
        
# Standalone code

if __name__ == '__main__':
    
    import time
    
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
    # ] # |B| = 15, nL = 5
    # B = [
    #     0b10011,
    #     0b01100,
    #     0b11000,
    #     0b00011,
    #     0b01001,
    #     0b10100,
    #     0b00110,
    #     0b01110
    # ] # |B| = 8, nL = 5
    # B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011] # Example from the article
    # B = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100, 0b0010, 0b0011] # 8-orbit
    # B = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100, 0b0010]
    # B = [6, 3, 1, 5, 0, 4, 2] # cost = 0
    B = [6, 2, 1, 0, 5]
    # B = [6,5]

    # B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011, 0b0000, 0b1111, 0b1011, 
    #      0b1101, 0b0110, 0b0010, 0b0101, 0b1000, 0b0001, 0b0111] # PROBLEM
    
    print(f"\nB = {[f'{b:0{len(bin(max(B)))-2}b}' for b in B]}")
    
    lxmixer = LXMixer(B, 3)
    # lxmixer = LXMixer(B, 4)
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
    
    print("\nnode_connectors =")
    for k, v in lxmixer.node_connectors.items():
        print(f"{k}")
        for neighbor, X in v.items():
            print(f"  <-> {neighbor} {X:0{lxmixer.nL}b}")

    print("\nOrbits (without projectors and costs):")
    for nodes, orbit in lxmixer.orbits.items():
        print(f"{nodes} : [{', '.join(f'{X:0{lxmixer.nL}b}' for X in orbit.Xs)}]")
    
    # """
    
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
        print(f"X costs = [{', '.join(str(X_cost) for X_cost in orbit.X_costs)}], Z cost = {orbit.Z_cost}")
    
    print("\nFinding best mixer...")
    start_time = time.time()
    lxmixer.find_best_mixer()
    best_combinations, best_Xs, best_Zs, best_cost = lxmixer.best_combinations, lxmixer.best_Xs, lxmixer.best_Zs, lxmixer.best_cost
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")
    
    print(f"\nFound {len(best_Xs)} best combinations of orbits with cost {best_cost}.")
    print("\nBest mixer:")
    print(f"[{', '.join(f'[{", ".join(f"[{", ".join(f"{x:0{lxmixer.nL}b}" for x in sub_Xs)}]" for sub_Xs in Xs)}]' for Xs in best_Xs)}]")
    print("\nBest projectors:")
    print(f"[{', '.join(f'[{", ".join(f"[{", ".join(f'{"+" if z[0] > 0 else "-"}{z[1]:0{lxmixer.nL}b}' for z in sub_Zs)}]" for sub_Zs in Zs)}]' for Zs in best_Zs)}]")    
   
    # plot_mixer_graph(lxmixer)
    # """