# import networkx as nx
# import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass, field
from itertools import combinations
import math
from tqdm import tqdm
import sys

from Stabilizer import *
from utils import ncnot, is_connected, split_into_suborbits, is_power_of_two, find_best_cost
from plot_mixers import draw_best_graphs, draw_mixer_graph

# TODO: Implement method for directed graphs (digraph=True). Only for visual representation.
# TODO: Implement method for reduced graphs (reduced=True)
# TODO: Implement blacklist and whitelist methods

@dataclass
class Orbit:
    """
    Class to store orbits and their properties.
    
    Attributes:
        Xs (Set[int]): Logical X operators.
        Zs (List[tuple[int, int]]): Projectors (Z operators).
        cost (int): Total cost (number of CNOTs required) of the orbit.
    """
    Xs: Set[int] = field(default_factory=list)
    Zs: List[Tuple[int, int]] = field(default_factory=list)
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
                
        processed_nodes = set()
        
        for seed in range(self.nB):
            #print(f"Processing seed {seed}")
            if seed in processed_nodes:
                #print(f"Seed {seed} already processed, skipping.")
                continue
            
            seed_Xs = list(self.node_connectors[seed].values())
            stack = []
            stack.append(([], seed_Xs, tuple([seed])))
            
            while stack:
                current_path, available_Xs, current_nodes = stack.pop()
                current_nodes_tuple = tuple(sorted(current_nodes))
                if len(current_nodes) == 2:
                    self.orbits[current_nodes_tuple] = Orbit(Xs=current_path)
                elif self.orbits.keys():
                    #print(self.orbits.keys())
                    for nodes in list(self.orbits.keys()):
                        if set(nodes).issubset(set(current_nodes)):
                            self.orbits[tuple(sorted(current_nodes))] = Orbit(Xs=current_path)
                            self.orbits.pop(nodes)
                            processed_nodes.update(current_nodes)
                    if not any(set(current_nodes).issubset(set(nodes)) for nodes in list(self.orbits.keys())):
                        self.orbits[current_nodes_tuple] = Orbit(Xs=current_path)
                        processed_nodes.update(current_nodes)
            
                if len(current_path) == len(seed_Xs):
                    continue
                
                for x, X in enumerate(available_Xs):
                    new_path = current_path + [X]
                    new_available = available_Xs[:x] + available_Xs[x+1:]
                    new_nodes = set()
                    for node in current_nodes:
                        new_nodes.update(n for edge in self.family_of_valid_graphs[X] for n in edge if node in edge)
                    
                    if len(new_nodes) == len(current_nodes) or not is_power_of_two(len(new_nodes)):
                        continue
                        
                    stack.append((new_path, new_available, tuple(sorted(new_nodes))))
                    
    def compute_costs(self):
        """
        Computes and updates the costs in the Orbit objects in `orbits`.
        """
        for nodes, orbit in self.orbits.items():
            best_Xs, best_cost = find_best_cost(orbit.Xs, orbit.Zs)
            orbit.Xs = best_Xs
            orbit.cost = best_cost
        
    def find_best_mixer(self):
        """ 
        Finds the best mixer based on the computed orbits, edges, minimal generating sets, projectors and costs.
        """
        
        if len(self.orbits.keys()) == 1:
            self.best_Xs = [[list(self.orbits.values())[0].Xs]]
            self.best_Zs = [[[(1,0)]]] #TODO
            self.best_cost = list(self.orbits.values())[0].cost
            self.best_combinations = [tuple(self.orbits.keys())]
            return
        
        N = range(2, len(self.orbits.keys())+1)
        for n in N:
            for combination in combinations(self.orbits.keys(), n):
                # time.sleep(0.05)
                if len(set([node for nodes in combination for node in nodes])) != self.nB:
                    #print(f"Combination {combination} does not cover all nodes in B, skipping.")
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
    B = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100, 0b0010, 0b0011] # 8-orbit
    # B = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100, 0b0010]
    # B = [6, 3, 1, 5, 0, 4, 2]
    # B = [6, 2, 1, 0, 5]
    # B = [6,5]

    # B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011, 0b0000, 0b1111, 0b1011, 
    #      0b1101, 0b0110, 0b0010, 0b0101, 0b1000, 0b0001, 0b0111] # PROBLEM
    
    print(f"\nB = {[f'{b:0{len(bin(max(B)))-2}b}' for b in B]}")
    
    # lxmixer = LXMixer(B, 3)
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
    for orbits in S.orbit_dictionary.values():
        print("Minimal generating sets:", orbits.Zs)
    end_time = time.time()
    print(f"\nTime: {end_time - start_time:.4f} s")
    print("______________")
    
    print("\nComputing projectors...")
    start_time = time.time()
    S.compute_projector_stabilizers()
    for orbits in S.orbit_dictionary.values():
        print("Projectors:", orbits.Zs)
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
   
    # """
    # """
   
    # draw_best_graphs(lxmixer)
    # fig, ax = plt.subplots()
    # draw_mixer_graph(ax, [list(lxmixer.orbits.keys())[0]], [list(lxmixer.orbits.values())[0].Xs], lxmixer, -0.1, r=0.1)
    # plt.show()
    # """