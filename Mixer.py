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
    def __init__(self, B, digraph=False, reduced=True, sort=False, blacklist=[], whitelist=None):
        self.setB(B, sort)
        
        self.digraph=digraph
        self.reduced=reduced
        self.base_solution_reduced = []
        self.base_solution = []
        self.base_cost=0
        self.base_nedges=0
        
        self.family_of_valid_graphs : Dict[int, List[Tuple[int,...]]] = {}
        # self.family_of_valid_graphs_flattened : Dict[int, List[int]]
        self.node_connectors = Dict[int, Dict[int, int]] # Maps nodes to connected nodes to the logical X operators that connect them
        
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

    def setB(self, B, sort:bool):
        """
        Sets the feasible set B for the mixer.

        Args:
            B (array-like[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        """
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
            B (list[int]): Feasible set of bitstrings (binary int representations) from the computational basis.
        
        Returns:
            Dict[int, List[Tuple[int,...]]]: A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator.
        """
        used_X = set()
        
        for i in range(self.nB):
            for j in range(i+1, self.nB):
                X_ij = self.B[i] ^ self.B[j]
                if X_ij not in used_X:
                    connected_indices = [k for k in range(self.nB) if (self.B[k] ^ X_ij) in self.B]
                    # self.family_of_valid_graphs_flattened[X_ij] = connected_indices
                    edges = set()
                    for k in connected_indices:
                        paired_val = self.B[k] ^ X_ij
                        paired_index = self.B.index(paired_val)
                        if paired_index > k:
                            edges.add((k, paired_index))
                    self.family_of_valid_graphs[X_ij] = edges
                    used_X.add(X_ij)
        
        # 2nd TRY: Using a dictionary to store edges with node values instead of a list
        # for i in range(self.nB):
        #     for j in range(i+1, self.nB):
        #         X_ij = self.B[i] ^ self.B[j]  # Logical X operator connecting nodes i and j
        #         if X_ij not in used_X:
        #             unconnected_mask = abs(1 - np.isin(self.B, X_ij^self.B))
        #             b1 = np.delete(self.B, unconnected_mask)
        #             b2 = b1 ^ X_ij
        #             self.family_of_valid_graphs[X_ij].append(frozenset(tuple(sorted(u, v)) for u, v in zip(b1,b2)))
        #             used_X.add(X_ij)
        
        # 1st TRY: Using a list to XGraph structs
        # graphs = []
        
        # for i in range(self.nB):
        #     for j in range(i+1, self.nB):
        #         edges = set()
        #         X_ij = self.B[i] ^ self.B[j]  # Logical X operator connecting nodes i and j
        #         if X_ij in used_X:
        #             continue
        #         unconnected_mask = abs(1 - np.isin(self.B, X_ij^self.B)) # Nodes in B that are not connected by X_ij
        #         b1 = np.delete(self.B, unconnected_mask)
        #         b2 = b1^X_ij
        #         # edges = list(zip(b1, b2))
        #         edges.add(tuple(sorted(u,v) for u, v in zip(b1, b2)))
        #         graphs.append(XGraph(X=X_ij, E = np.array(edges))) # G can be a list/array
                
        #         used_X.add(X_ij)
        
        # self.family_of_valid_graphs = graphs
        
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
        failed_prefixes = set() # Set to store prefixes of failed paths to avoid using these combinations again

        for seed in range(self.nB):
            if seed in processed_nodes: # Check if the seed node has already been processed
                continue
            
            seed_Xs = list(self.node_connectors[seed].values()) # All X operators connected to the seed node
            
            stack = []
            initial_available = set(seed_Xs)
            stack.append(([], initial_available, set([seed]))) # Initializes depth-first search with an empty path, available X operators and the seed node
            
            while stack:
                current_path, available_Xs, current_nodes = stack.pop() # Processes each state in last-in-first-out order
                # current_path: sequence of X operators applied so far
                # available_Xs: set of X operators that can still be applied
                # current_nodes: current nodes in orbit
                
                path_tuple = tuple(sorted(current_path))
                if path_tuple in failed_prefixes:
                    continue # Skip paths if combinations of operators have already failed
                
                if len(current_nodes) >= 2:
                    orbit_tuple = tuple(sorted(current_path))
                    if not any(orbit_tuple == tuple(sorted(existing)) for existing in self.orbits): # Check if the orbit is already recorded
                        self.orbits.append(list(current_path))
                        self.nodes.append(list(current_nodes))
                        
                if len(current_path) == len(seed_Xs): # If all available X operators have been used, go to a new branch
                    continue
                
                for X in available_Xs:
                    new_path = current_path + [X]
                    new_nodes = list({node for u, v in self.family_of_valid_graphs[X] if u in current_nodes and v in current_nodes for node in (u, v)})
                    
                    if not new_nodes: 
                        failed_prefixes.add(tuple(sorted(new_path)))
                        continue
                    
                    # Prepare next step
                    new_available = set(available_Xs) - {X}
                    stack.append((new_path, new_available, new_nodes))

        # 3rd TRY ###########################
            # for X_sequence in permutations(seed_Xs, self.nB-1):
            #     nodes = self.family_of_valid_graphs_flattened[X_sequence[0]]
            #     orbit = set()
                
            #     if len(nodes) == 2 : # Check if the orbit is trivial
            #         self.orbits.append([X_sequence[0]])
            #         self.nodes.append[nodes]
            #         continue
                
            #     sequence_length = 0
            #     for X in X_sequence:
            #         new_nodes = [u in nodes and v in nodes for u,v in self.family_of_valid_graphs[X]]
            #         if len(new_nodes) == 0:
            #             break
            #         nodes = new_nodes
            #         orbit.update(X)
            #         sequence_length += 1
                
            #     self.orbits.append(list(orbit))
            #     self.nodes.append(nodes)

        # 2nd TRY ###########################
        # queue = [current_node]
        
        # while len(processed_nodes) < self.nB:
        #     if current_node in processed_nodes:
        #         current_node = min(set(range(self.nB)) - processed_nodes)
        #         continue
        
        #     orbit_Xs = set()
        #     orbit_Vs = {current_node}
        #     orbit_Es = set()
            
        #     queue = [current_node]
            
        #     while queue:
        #         node = queue.pop(0)
                
        #         for neighbor, X in self.node_connectors[node] :
        #             if neighbor not in orbit_Vs:
        #                 orbit_Xs.add(X)
        #                 orbit_Vs.add(neighbor)
        #                 orbit_Es.add(tuple(sorted((node, neighbor))))
                        
        #                 queue.append(neighbor)
                
        #     self.orbits.append(list(orbit_Xs))
        #     self.nodes.append(list(orbit_Vs))
            
        #     processed_nodes.update(orbit_Vs)
            
        #     if len(processed_nodes) < self.nB:
        #         current_node = min(set(range(self.nB)) - processed_nodes)
            
        # 1st TRY ###########################
        # iteration = 0
        # included_nodes = set([0])

        # while included_nodes != set(self.B):
        #     Xs = self.node_connectors[current_node][iteration:]
        #     num_X = len(Xs)
        #     Es = [self.family_of_valid_graphs[X] for X in Xs]
        #     for i, E in enumerate(Es):
        #         V = [i for t in E for i in t]
        #         if len(E) == 1:
        #             self.orbits.append([Xs[i]])
        #             self.edges.append(E)
        #             continue
        #         branch_length = 1
        #         orbit_Xs = [Xs[0]]
        #         current_X = orbit_Xs[-1]
        #         current_edges = E
        #         while len(current_edges) > 1 and branch_length < num_X:
        #             for e, edge in enumerate(current_edges):
        #                 if edge[0] ^ current_X not in V and edge[1] ^ current_X not in V:
        #                     del current_edges[e]
                    
        #             orbit_Xs.append(Xs((i+branch_length) % num_X))
        #             current_X ^= orbit_Xs[-1]
                
        #         self.edges.append(current_edges)
        #         self.orbits.append([orbit_Xs])
        #         included_nodes.update([i for t in current_edges for i in t])
                
        # iteration += 1
        # current_node = (i for i in included_nodes not in set(self.B))[0]
                    
            
                

    # def compute_all_subgraphs(self, xgraph_list:list[XGraph]):
    #     """
    #     Computes all subgraphs of a given list of XGraph.
        
    #     Args:
    #         graph_list (list[XGraphs]]): List of XGraph structs to divide into subgraphs.
    #     """
    #     for xG in xgraph_list:
    #         edges = list(xG.G.edges)
    #         for r in range(1, len(edges) + 1):
    #             if is_power_of_two(r):
    #                 for subset in combinations(edges, r):
    #                     subG = nx.Graph()
    #                     subG.add_edges_from(subset)
    #                     if check_if_orbit(subG.nodes)[0]:
    #                         xG.subgraphs.append(subG)
                        
    #                         xG.projectors.append(np.array(compute_restricted_projector_stabilizer(subG.nodes(), self.nL)))
    #                         xG.subgraph_costs.append(np.sum(ncnot(xG.projectors[-1])))
                        
    #     # pass-by-reference  
    
    # def check_all_combinations(self, graph_list:list[XGraph], N=None, combine_subgraphs=True):
    #     """
    #     Checks all combinations of (sub)graphs in the given list of XGraph.
    #     If n is specified, only combinations of size n are checked.

    #     Args:
    #         graph_list (list[XGraph]): List of XGraph structs to check combinations from.
    #         N (int, optional): Size of combinations. Defaults to None.

    #     Returns:
    #         list[XGraph]: List of XGraph structs that are connected.
    #     """
    #     best_mixer = None
        
    #     if N is None: N = range(1,len(graph_list)+1)
    #     else: N = range(1,N+1)
        
    #     options = []
    #     for xG in graph_list:
    #         if combine_subgraphs:
    #             options.append([(xG.X, subgraph, projectors_i, subgraph_cost) for subgraph, projectors_i, subgraph_cost in zip(xG.subgraphs, xG.projectors, xG.subgraph_costs)] + [(None, nx.Graph(), None, None)])
    #         else: options.append([(xG.X, xG.G, None, xG.cost), (None, nx.Graph(), None, None)])
        
    #     current_best_cost = float('inf')
        
    #     for n in N:
    #         for subset in combinations(options,n): # Unordered combinations of size n
    #             for combination in product(*subset):
    #                 combined_G = nx.Graph()
    #                 Xs = []
    #                 projectors = []
    #                 combined_cost = 0
    #                 for X_i, subG, projectors_i, cost in combination:
    #                     combined_G = nx.compose(combined_G, subG)
    #                     if X_i is not None: Xs.append(X_i)
    #                     projectors.append(projectors_i)
    #                     combined_cost += cost
    #                 #  Replace if lower cost found               
    #                 if combined_G.number_of_nodes == len(self.B) and combined_G.number_of_edges == (len(self.B)-1): # Easy check
    #                     if combined_cost < current_best_cost:
    #                         if nx.is_connected(combined_G): # Actual check
    #                             best_mixer = CombinedXGraph(
    #                                     Xs=Xs, 
    #                                     G=combined_G, 
    #                                     projectors=projectors,
    #                                     cost = combined_cost
    #                                     )
    #                             current_best_cost = combined_cost
    #     return best_mixer

