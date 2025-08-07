import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass, field
from itertools import combinations, product
import math
from tqdm import tqdm

from Stabilizer import *
from utils import is_connected, is_power_of_two, find_best_cost
from plotting.plot_mixers import draw_best_graphs, draw_mixer_graph, draw_family_of_valid_graphs, X0

# TODO: Finish implementing the "semi_restricted_suborbits" method. The current implementation runs fine, and even yields a lower cost than the
# "all_suborbits" method (as expected), but wrongly excludes other nodes in the main orbit that are connected by the same operator as a suborbit
# (and thus probably finds a higher than necessary cost solution).

@dataclass
class Suborbit:
    Xs: List[int] = field(default_factory=list)
    cost: int = float('inf')

@dataclass
class Orbit:
    """
    Class to store orbits and their properties.
    
    Attributes:
        Xs (List[int]): Logical X operators.
        Zs (List[Tuple[int, int]]): Projectors (Z operators). First element in the tuple is the signe (+1/-1) and the second element the operator (int representation).
        cost (int): Total cost (number of CNOTs required) of the orbit.
    """
    Xs: List[int] = field(default_factory=list)
    Zs: List[Tuple[int, int]] = field(default_factory=list)
    cost: int = float('inf') # Total cost of the orbit, initialized to infinity.
    
    suborbits: Dict[Tuple[int,...], Suborbit] = field(default_factory=dict)

class LXMixer:
    """
    Logical X mixer for the QAOA.
    
    The mixer is based on logical X operators connecting nodes in a feasible set B, the span of which is the feasible solution space
    of the QAOA problem. The mixer computes the family of valid graphs, orbit subgraphs within this family and finds the minimal cost combinations of
    subgraphs (with their respective logical X operators and projectors) that connect all nodes in the feasible set B.
    
    Attributes:
        B (List[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        nB (int): Number of elements in the feasible set B.
        nL (int): Number of qubits. 
        sort (bool): Whether to sort the feasible set B (default: False).
        method (str): Method to use for finding the best mixer, either "largest_orbits", "all_suborbits" or "semi_restricted_suborbits" (default: "largest_orbits").
        
        family_of_valid_graphs (Dict[int, List[Tuple[int,...]]]): A dictionary mapping logical X operators (int) to edges (tuples of node indices) connected by the operator X_ij : (i,j).
        node_connectors (Dict[Dict[int, int]]): Maps for each node, another node from B and the logical X operator that connects them i : j : X_ij.
        orbits (Dict[Tuple[int,...], Orbit]): A dictionary mapping node tuples to Orbit objects, containing the logical X operators, the Z operators representing the orbit and their total cost.
        
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
    def __init__(self, B, nL, sort=False, method="largest_orbits"):
        """
        Initializes the LXMixer with the feasible set B and number of logical qubits nL.

        Args:
            B (List[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
            nL (int): Number of quibits.
            sort (bool): Whether to sort the feasible set B (default: False).
        """
        self.setB(B, nL, sort)
        
        if method not in ["largest_orbits", "all_suborbits", "semi_restricted_suborbits"]:
            raise ValueError("Method must be one of 'largest_orbits', 'all_suborbits' or 'semi_restricted_suborbits'.")
        self.method = method
        
        # Prepare empty attributes
        self.family_of_valid_graphs : Dict[int, List[Tuple[int,...]]] = {} # Maps logical X operators to edges (tuples of node indices) connected by the operator X_ij : (i,j).
        self.node_connectors : Dict[int, Dict[int, int]] = {} # Maps for each node, another node from B and the logical X operator that connects them i : j : X_ij.
        for i in range(self.nB): self.node_connectors[i] = {}
        self.orbits : Dict[Tuple[int,...], Orbit] = {}
        self.best_Xs = []
        self.best_Zs = []
        self.best_combinations = []
        self.best_cost: int = float('inf') # Set best cost to infinity initially.

    def setB(self, B, nL, sort:bool):
        """
        Sets the feasible set B for the mixer.

        Args:
            B (iterable[int]): Feasible set of bitstrings (int representations) from the computational basis.
        """
        if isinstance(B, set):
            B = list(B)
        elif isinstance(B, list):
            seen = set()
            B = [x for x in B if not (x in seen or seen.add(x))] # Remove duplicates while preserving order
        else:
            raise TypeError("B must be a list or a set.")
        self.nB = len(B) # |B|
        if sort: # Sort B in ascending order if sort is True
            B=sorted(B, key=lambda x: int(x, 2))
            
        if len(B) < 2:
            raise Exception("B must contain at least two elements.")

        self.nL = nL
        for b in B:
            if b >= (1 << self.nL):
                raise Exception(f"Entry {b} exceeds {self.nL} bits.")

        self.B = B
            
    def compute_family_of_valid_graphs(self):
        """
        Computes the family of valid graphs for the feasible set B.
        """        
        for i in range(self.nB):
            for j in range(i+1, self.nB):
                X_ij = self.B[i] ^ self.B[j] # XOR operation to find the logical X operator connecting nodes i and j.
                if X_ij not in self.family_of_valid_graphs.keys():
                    self.family_of_valid_graphs[X_ij] = [(i,j)]
                else: 
                    self.family_of_valid_graphs[X_ij].append((i,j))
        
    def compute_all_orbits(self):
        """
        Computes all group-generated closed sets ("orbits") in the family of valid graphs using a depth-first search algorithm.
        Results stored in the `orbits` attribute as a dictionary mapping node tuples to Orbit objects, updating their logical X operators `Xs`.
        """    
                
        # Maps for each node the |B|-1 logical X operators that connects it to the other nodes i : j : X_ij.
        # Is used for extracting all the logical X operators that connect to each seed node (`seed_Xs =` ...).
        for X, E in self.family_of_valid_graphs.items():
            for i, j in E:
                self.node_connectors[i][j] = X
                self.node_connectors[j][i] = X
                
        processed_nodes = set() # Set to store nodes that have already been associated with a non-trivial orbit.
        
        for seed in range(self.nB): # Iterate over each node in B as a seed node.
            if seed in processed_nodes and self.method == "largest_orbits": # Skip if the seed node is already processed (unless wanting all suborbits).
                continue
            
            seed_Xs = list(self.node_connectors[seed].values()) # The |B|-1 logical X operators that connect to the seed node.
            stack = []
            stack.append(([], seed_Xs, tuple([seed]))) # Append the first path to the stack: empty path, available X operators and the seed node as a tuple.
            
            while stack:
                current_path, available_Xs, current_nodes = stack.pop() # Backtracking step. Retrieve the last path from which we arrived at the last "leaf".
                """
                current_path: List of logical X operators that form the current orbit.
                available_Xs: List of logical X operators that can be added to the current path, i.e. those that have not been used yet.
                current_nodes: Tuple of nodes that are connected by the current path.
                """
                if len(current_nodes) > 1:
                    current_nodes_tuple = tuple(sorted(current_nodes)) # Sort the nodes for consistent key representation.
                    
                    subset_found = any(set(current_nodes).issubset(set(nodes)) for nodes in list(self.orbits.keys())) # If the current nodes are a subset of any existing orbit.
                    
                    # If orbits doesn't contain anything, or if the current nodes are not a subset of any existing orbit, or if we want to find all suborbits, add the current nodes as a new orbit.
                    if not self.orbits.keys() or not subset_found or self.method == "all_suborbits":
                        self.orbits[current_nodes_tuple] = Orbit(Xs=current_path)
                        
                        # If we only want the largest orbits in the main dictionary, remove sets that are subsets of the new orbit.
                        if self.method != "all_suborbits":
                            for nodes in list(self.orbits.keys()):
                                if set(nodes) < set(current_nodes):
                                    if self.method == "semi_restricted_suborbits": # If we want to keep suborbit information (semi-restricted suborbits) add suborbits to their respective main orbits.
                                        self.orbits[current_nodes_tuple].suborbits[nodes] = Suborbit(Xs=self.orbits[nodes].Xs)
                                    self.orbits.pop(nodes)
                                    
                        if len(current_nodes) > 2: processed_nodes.update(current_nodes) # Mark nodes in non-trivial orbits as processed. Will be skipped when these seeds are reached.
                
                    elif self.method == "semi_restricted_suborbits": # If we want to keep suborbit information (semi-restricted suborbits), and some orbit contains the current nodes, add the current nodes as a suborbit to this larger orbit.
                        for nodes in self.orbits.keys():
                            if set(nodes) > set(current_nodes):
                                self.orbits[nodes].suborbits[current_nodes_tuple] = Suborbit(Xs=current_path)
                
                for x, X in enumerate(available_Xs): # Iterate over the available X operators to choose from at this point in the tree. 
                    new_path = current_path + [X] # Add the new X operator to the path.
                    new_available = available_Xs[:x] + available_Xs[x+1:] # Remove the new X operator from the available operators.
                    new_nodes = set() # Set of nodes that are connected to the current nodes by the new X operator.
                    
                    for node in current_nodes: 
                        new_nodes.update(n for edge in self.family_of_valid_graphs[X] for n in edge if node in edge) # Add new nodes that are connected to the current nodes by the new X operator.
                    
                    if int(math.log2(len(new_nodes))) != len(new_path) or not is_power_of_two(len(new_nodes)): # If the new X operator doesn't connect nodes that are a power of two or doesn't increase the numeber of nodes, don't add it to the stack.
                        continue
                        
                    stack.append((new_path, new_available, tuple(sorted(new_nodes)))) # Add a valid path to the stack.
        
        if self.method == "semi_restricted_suborbits":
            for nodes, orbit in self.orbits.items():
                orbit.suborbits[nodes] = Suborbit(Xs=orbit.Xs) # Add the orbit itself as a suborbit to itself (simpler when considering combinations).
                    
    def compute_costs(self):
        """
        Computes and updates the costs in the Orbit objects in `orbits`.
        
        Iterates through each orbit and chooses the combination of logical X operators that yield the lowest cost with the Z operators.
        """        
        for nodes, orbit in self.orbits.items():
            best_Xs, best_cost = find_best_cost(orbit.Xs, orbit.Zs) # Finds the best combination of log2(# of nodes) logical X operators that yield the lowest cost with the projectors.
            orbit.Xs = best_Xs
            orbit.cost = best_cost
            for subnodes, suborbit in orbit.suborbits.items():
                if subnodes == nodes: 
                    orbit.suborbits[subnodes].Xs = orbit.Xs
                    orbit.suborbits[subnodes].cost = orbit.cost
                    continue
                best_Xs, best_cost = find_best_cost(suborbit.Xs, orbit.Zs)
                suborbit.Xs = best_Xs
                suborbit.cost = best_cost
        
    def find_best_mixer(self):
        """ 
        Finds the best mixer and updates the `best_Xs`, `best_Zs`, `best_combinations` and `best_cost` attributes.
        """ 
        # If all the nodes in B are connected by a single orbit, this is the best mixer.
        B_orbit = tuple(range(self.nB)) # All node indeces in B.
        if B_orbit in self.orbits.keys():
            self.best_Xs = [[self.orbits[B_orbit].Xs]]
            self.best_Zs = [[self.orbits[B_orbit].Zs]]
            self.best_cost = self.orbits[B_orbit].cost
            self.best_combinations = [[B_orbit]]
            return
        
        best_main_combinations = [] # List to store the best main combinations of orbits (if using semi-restricted suborbits).
        N = range(2, min([self.nB, len(self.orbits)+1])) # Range of the number of orbits to combine. Goes from 2 to |B|-1 (worst-case is a chain).
        for n in N:
            if self.method == "semi_restricted_suborbits": # If using semi-restricted suborbits, we need to consider combinations of main orbits and their suborbits.
                for main_combination in tqdm(combinations(self.orbits.keys(), n), desc=f"Main combo size {n}/{self.nB-1}"):
                    if len({node for nodes in main_combination for node in nodes}) != self.nB or not is_connected(main_combination): # If the combination doesn't cover all nodes in B or is unconnected, skip.
                        continue
                    for combination in product(*[self.orbits[main_nodes].suborbits.keys() for main_nodes in main_combination]): # Cartesian product of the main orbits' suborbits.
                        if len({node for nodes in combination for node in nodes}) != self.nB or not is_connected(combination): # If the combination doesn't cover all nodes in B or is unconnected, skip.
                            continue
                        cost = 0
                        for main_nodes, sub_nodes in zip(main_combination, combination):
                            cost += self.orbits[main_nodes].suborbits[sub_nodes].cost
                            if cost > self.best_cost: # If the cost at any point exceeds the best cost, go to the next combination.
                                break
                        if cost < self.best_cost: # If the cost is strictly lower that the best cost, update the best cost and combination.
                            self.best_cost = cost
                            best_main_combinations = [main_combination]
                            self.best_combinations = [combination]
                        elif cost == self.best_cost: # If the cost is equal to the best cost, add the combination to the list of best combinations.
                            best_main_combinations.append(main_combination)
                            self.best_combinations.append(combination)
            else:
                for combination in tqdm(combinations(self.orbits.keys(), n), desc=f"Combo size {n}/{self.nB-1}"):
                    if len({node for nodes in combination for node in nodes}) != self.nB: # If the combination does not cover all nodes in B, skip it.
                        continue
                    if not is_connected(combination): # If the combination doesn't connect all nodes or is unconnected, skip.
                        continue
                    cost = 0
                    for orbit_nodes in combination:
                        cost += self.orbits[orbit_nodes].cost # Add the cost each orbit in the combination.
                        if cost > self.best_cost: # If the cost at any point exceeds the best cost, go to the next combination.
                            break
                    if cost < self.best_cost: # If the cost is strictly lower that the best cost, update the best cost and combination.
                        self.best_cost = cost
                        self.best_combinations = [combination]
                    elif cost == self.best_cost: # If the cost is equal to the best cost, add the combination to the list of best combinations.
                        self.best_combinations.append(combination)
        
        # Add the corresponding X operators and projectors to the best combination(s).
        if self.method == "semi_restricted_suborbits":
            self.best_Xs = [
                [self.orbits[main_nodes].suborbits[sub_nodes].Xs for main_nodes, sub_nodes in zip(main_comb, sub_comb)]
                for main_comb, sub_comb in zip(best_main_combinations, self.best_combinations)
            ]
            self.best_Zs = [
                [self.orbits[main_nodes].Zs for main_nodes in main_comb]
                for main_comb in best_main_combinations
            ]
        else:
            self.best_Xs = [[self.orbits[orbit_nodes].Xs for orbit_nodes in combination] for combination in self.best_combinations]
            self.best_Zs = [[self.orbits[orbit_nodes].Zs for orbit_nodes in combination] for combination in self.best_combinations]
        return
        
# Standalone code, e.g. example usage and testing.
if __name__ == '__main__':
    import time # For measuring execution time.
    
    # nL = 3
    # B = [6,5] # nB = 2
    # B = [6, 2, 1, 0, 5] # nB = 5
    # B = [6, 3, 1, 5, 0, 4, 2] # nB = 7
    # B = [0, 1, 2, 3, 4, 5, 6, 7] # bB = 8, whole space, 8-orbit
    
    # nL = 4
    # B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011] # nB = 5, example from the article
    # B = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100] # nB = 6
    # B = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100, 0b0010] # nB = 7
    B = [0b0000, 0b1111, 0b0001, 0b1101, 0b1110, 0b1100, 0b0010, 0b0011] # nB = 8, 8-orbit
    # B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011, 0b0000, 0b1111, 0b1011, 
    #      0b1101, 0b0110, 0b0010, 0b0101, 0b1000, 0b0001, 0b0111] # nB = 15
    B = [0,1,2,3,4,5,6,7]
    # nL = 5
    # B = [0b10011, 0b01100, 0b11000, 0b00011,
    #     0b01001, 0b10100, 0b00110, 0b01110] # nB = 8
    # B = [0b00001, 0b00010, 0b00100,
    #     0b01000, 0b10000, 0b00011,
    #     0b00101, 0b00110, 0b01001,
    #     0b01010, 0b01100, 0b10001, 
    #     0b10010, 0b10100, 0b11000] # nB = 15
    
    print(f"\nB = {[f"{b:0{len(bin(max(B)))-2}b}" for b in B]}") # Print B in binary format.
    
    # Initialize the LXMixer with the feasible set B and number of logical qubits nL.
    # lxmixer = LXMixer(B, 3)
    lxmixer = LXMixer(B, 4, method="semi_restricted_suborbits")
    # lxmixer = LXMixer(B, 4, method="all_suborbits")
    # lxmixer = LXMixer(B, 4, method="largest_orbits")
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
        print(f"{nodes} : Xs = [{', '.join(f'{X:0{lxmixer.nL}b}' for X in orbit.Xs)}]")
        for suborbit_nodes, suborbit in orbit.suborbits.items():
            print(f"  {suborbit_nodes} : Xs = [{', '.join(f'{X:0{lxmixer.nL}b}' for X in suborbit.Xs)}]")
    # """
    
    S = Stabilizer(lxmixer.B, lxmixer.nL, lxmixer.orbits) # Initialize the Stabilizer object with the feasible set B, number of logical qubits nL and orbits dictionary.
    
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
        print(f"{nodes} : Xs = [{', '.join(f'{X:0{lxmixer.nL}b}' for X in orbit.Xs)}], Zs = [{', '.join(f'{"+" if Z[0] == 1 else "-"}{Z[1]:0{lxmixer.nL}b}' for Z in orbit.Zs if len(Z) == 2)}], cost = {orbit.cost}")
        for suborbit_nodes, suborbit in orbit.suborbits.items():
            print(f"  {suborbit_nodes} : Xs = [{', '.join(f'{X:0{lxmixer.nL}b}' for X in suborbit.Xs)}], cost = {suborbit.cost}")
    
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
    
    # Draw family of valid graphs.
    # draw_family_of_valid_graphs(lxmixer, lw=1.5, r=0.2, group_size=2)
   
    # Draw the best mixer graph.
    draw_best_graphs(lxmixer, r=0.15, lw =1.25)
    
    # Draw specified orbit(s).
    # fig, ax = plt.subplots()
    # draw_mixer_graph(ax, [list(lxmixer.orbits.keys())[0]], [list(lxmixer.orbits.values())[0].Xs], lxmixer, x=X0, r=0.1)
    # """
    
    plt.show()
    
    # If wanting to save plots, use the `saveas` parameter in the drawing functions.