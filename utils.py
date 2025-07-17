import math
import networkx as nx
from itertools import combinations
from functools import reduce
import operator

is_power_of_two = lambda x: (x > 0) and (x & (x - 1)) == 0

def ncnot(P) :
    """
    Calculate the number of CNOT gates required to implement a Pauli string.

    Args:
        P (int): Pauli string represented as a binary integer, where each 1 represents a qubit that is acted upon by a Pauli operator (X or Z).

    Returns:
        int: Number of CNOT gates required to implement the Pauli string.
    """
    ncnot = P.bit_count()
    return (ncnot > 1)*2*(ncnot - 1)

def pauli_int_to_str(P, operator):
    P = str(P)
    P.replace("0", "I")
    if operator == "X":
        P = P.replace("1", "X")
    elif operator == "Z":
        P = P.replace("1", "Z")
    else:
        raise ValueError("Operator must be 'X' or 'Z'.")
    return P

def parity(n):
    #using Brian Kernighan's algorithm to check parity (commutation/anti-commutation), parity = 1 if even (commutes), and parity = -1 if odd (anti-commutes)
    parity = 0
    while n:
        parity ^= 1
        n &= n - 1 
    if parity == 0:
        parity = 1
    else:
        parity = -1
    return parity

def convert_to_binary_string(int_values, n):
    """
    Convert integers to binary strings with n bits.
    
    - If int_values is a single int: return binary string
    - If it's a list: return list of binary strings or (bin, bin) tuples
    - If it's a dict: return dict with binary string keys and values
    """
    
    if isinstance(int_values, int):
        return format(int_values, f'0{n}b')

    elif isinstance(int_values, list):
        result = []
        for i in int_values:
            if isinstance(i, int):
                result.append(format(i, f'0{n}b'))
            elif isinstance(i, tuple) and len(i) == 2:
                if all(isinstance(x, int) for x in i):
                    result.append((format(i[0], f'0{n}b'), format(i[1], f'0{n}b')))
        return result

    elif isinstance(int_values, dict):
        result = {}
        for key, val in int_values.items():
            if isinstance(key, int) and isinstance(val, int):
                result[format(key, f'0{n}b')] = format(val, f'0{n}b')
        return result

    else:
        raise TypeError("Input must be int, list, or dict of ints")

def is_connected(orbits):
    # Short circuit evaluation
    for orbit in orbits:
        if not any(set(orbit).intersection(set(other_orbit)) for other_orbit in orbits if other_orbit != orbit):
            return False
    return True

def split_into_suborbits(family_of_valid_graphs, operators, nodes = None):
    """
    Splits a set of nodes into suborbits based on the provided operators and family of valid graphs.

    Args:
        family_of_valid_graphs (Dict[int, List[Tuple[int,...]]]): A dictionary mapping logical X operators (int representations) to edges (tuples of node indices) connected by the operator. 
        operators (list): A list of logical X operators (int representations) that makes up the suborbits.

    Returns:
        List[set()]: A list of suborbits in sets with the nodes that are connected by the operators. 
    """
    suborbit_size = len(operators) + 1
    
    # Finds how many X operators are needed to cover the suborbit
    least_number_of_operators = int(math.log2(suborbit_size))

    all_edges = []

    # Iterates through needed number of operators
    for operator in operators[:least_number_of_operators]:
        list_of_edges = family_of_valid_graphs[operator] #TODO can these nodes in list of edges be outside the nodes for the orbit?
        # Only include edges where both nodes are in the specified nodes set
        all_edges += [edge for edge in list_of_edges if all(node in nodes for node in edge)]  
        #all_edges += list_of_edges TODO option 2: could it be more efficeent to filter out solutions later?
    
    G = nx.Graph()
    G.add_edges_from(all_edges)

    suborbits = list(nx.connected_components(G))
    #suborbits = [suborbit for suborbit in nx.connected_components(G) if len(suborbit) == suborbit_size] TODO option 2

    return suborbits

if __name__ == '__main__':
    family_of_valid_graphs_1 = {0b0010 : [(0, 1), (2, 7), (3, 9), (4, 13), (5, 10), (6, 8), (11, 14)],
    0b0111 : [(0, 2), (1, 7), (3, 4), (5, 14), (6, 12), (9, 13), (10, 11)],
    0b1010 : [(0, 3), (1, 9), (2, 4), (6, 11), (7, 13), (8, 14), (10, 12)],
    0b1101 : [(0, 4), (1, 13), (2, 3), (5, 8), (6, 10), (7, 9), (11, 12)],
    0b1110 : [(0, 5), (1, 10), (2, 14), (4, 8), (6, 13), (7, 11), (9, 12)],
    0b0001 : [(0, 6), (1, 8), (2, 12), (3, 11), (4, 10), (5, 13), (9, 14)],
    0b0101 : [(0, 7), (1, 2), (3, 13), (4, 9), (5, 11), (8, 12), (10, 14)],
    0b0011 : [(0, 8), (1, 6), (3, 14), (4, 5), (7, 12), (9, 11), (10, 13)],
    0b1000 : [(0, 9), (1, 3), (2, 13), (4, 7), (5, 12), (6, 14), (8, 11)],
    0b1100 : [(0, 10), (1, 5), (2, 11), (3, 12), (4, 6), (7, 14), (8, 13)],
    0b1011 : [(0, 11), (1, 14), (2, 10), (3, 6), (4, 12), (5, 7), (8, 9)],
    0b0110 : [(0, 12), (2, 6), (3, 10), (4, 11), (5, 9), (7, 8), (13, 14)],
    0b1111 : [(0, 13), (1, 4), (2, 9), (3, 7), (5, 6), (8, 10), (12, 14)],
    0b1001 : [(0, 14), (1, 11), (2, 5), (3, 8), (6, 9), (7, 10), (12, 13)],
    0b0100 : [(1, 12), (2, 8), (3, 5), (4, 14), (6, 7), (9, 10), (11, 13)]}
    operators_1 = [0b1001, 0b0110, 0b111]
    nodes_1 = (0, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 14)

    operators_2 = [0b0010, 0b0110, 0b1000, 0b1010, 0b1100, 0b1110]
    nodes_2 = (2, 6, 7, 8)
    
    suborbits = split_into_suborbits(family_of_valid_graphs=family_of_valid_graphs_1, operators=operators_2, nodes=nodes_2)
    #print("suborbits: ", suborbits)

def find_best_cost(Xs, Zs):#orbit):
    """
    takes in a list of Xs that are log2(nodes) and that generates an orbit. Also takes in a list of Zs that corresponds to the orbit.
    """
    # Xs = orbit.Xs 
    # Zs = orbit.Zs #remember to change from (1, string) to only string...
    all_x_operators = []
    n = len(Xs)
    
    # Generate all combinations of Xs and their corresponding hats
    for r in range(1, n + 1):  # Start from 1 to include single elements
        for combo in combinations(Xs, r):
            hat = reduce(operator.xor, combo)
            all_x_operators.append([combo, hat])
    
    all_costs = {}

    for used_Xs, X_combos in all_x_operators:
        total_cost = 0
        for Z in Zs:
            cost = ncnot(X_combos ^ Z)
            total_cost += cost
        
        all_costs[used_Xs] = total_cost
    
    """
    # TODO Check this part!!! might be able to do it more efficiently
    # The Xs we demand will be in the solution somehow (to make sure its an orbit)
    required_set = set(Xs)
    best = None
    min_cost = float('inf')
    all_costs_as_list = list(all_costs.items())
    for group in combinations(all_costs_as_list, n):
        covered_1 = set()
        total_cost = 0

        for combo, cost in group:
            covered_1.update(combo)
            total_cost += cost

        if covered_1 >= required_set and total_cost < min_cost:
            best = group
            min_cost = total_cost
    """

    # New version, hopefully more efficient
    best_Xs = []
    best_cost = 0
    covered = set()
    required = set(Xs)
    maybe_later = []
    
    while len(best_Xs) < n:
        lowest_cost = min(all_costs.values())
        keys = [k for k, v in all_costs.items() if v == lowest_cost]
        #print("this is the key:", keys)
        #print("this is the lowest cost:", lowest_cost)
        
        # iterate through the keys with the lowest cost
        for key in keys:
            # Checks that either the key adds to the subset or that it is already covered (i.e. that we are actually creating an orbit)
            if (not set(key).issubset(covered)) or (required == covered):  # If key has *any* uncovered elements
                # Checks that if it is already covered, we use the lowest cost from maybe_later
                if required == covered:
                    new_key_and_cost = maybe_later.pop(0) if maybe_later else [key, lowest_cost]
                    best_Xs.append(new_key_and_cost[0])
                    best_cost += new_key_and_cost[1]
                else:
                    covered.update(key)
                    best_Xs.append(key)
                    best_cost += lowest_cost
                
                if len(best_Xs) == n:
                    break
            
                del all_costs[key]
            else:
                del all_costs[key]  
                maybe_later.append([key, lowest_cost])


    return best_Xs, best_cost#, list(best), min_cost

if __name__ == '__main__':
    results = find_best_cost([0b0010, 0b0110, 0b1000, 0b1001, 0b1111, 0b0000], [0b0010, 0b0110, 0b1000, 0b1010, 0b1100, 0b1110])

    print("Best combo of Xs (heuristic):", results[0],"\nBest cost (heuristic):", results[1], "\nBest combo of Xs (exact):", results[2], "\nBest cost (exact):", results[3])