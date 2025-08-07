import math
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

def pauli_int_to_str(P, nL, operator="X"):
    """
    Converts a Pauli string represented as an integer to its string representation.

    Args:
        P (int): Pauli string.
        nL (int): Number of qubits.
        operator (str, optional): Whether the 1s represent "X" or "Z" operators. Defaults to "X".

    Raises:
        ValueError: If the operator is not "X" or "Z".

    Returns:
        str: Formatted Pauli string.
    """
    P = f"{P:0{nL}b}"
    P = P.replace("0", "I")
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
    - If int_values is a list: return list of binary strings or (bin, bin) tuples
    - If int_values is a dict: return dict with binary string keys and values
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
    """
    Checks if the provided orbits are connected, meaning that there is at least one node in common between any two orbits.

    Args:
        orbits (List[Tuple[int, ...]]): Node indices representing the group generated sets.

    Returns:
        bool: True if the orbits are connected, False otherwise.
    """
    # Short circuit evaluation
    if len(orbits) == 1: return True # Return True if there is only one orbit (connecting all nodes)
    for orbit in orbits:
        if not any(set(orbit).intersection(set(other_orbit)) for other_orbit in orbits if other_orbit != orbit):
            return False
    return True


def find_best_cost(Xs, Zs_operators):
    """
    Finds the best combination of logical X operators that generates an orbit, minimizing the cost function.
    
    Args:
        Xs (List[int]): A list of logical X operators (int representations) that generates an orbit.
        Zs_operators (List[Tuple[int, str]]): A list of tuples representing the Z-operators, where the first element is the sign(+1/-1) and the second is the int representation of the oeprator.
    
    Return:
        List[int]: A list of the best logical X operators (int representations).
    """
    all_x_operators = []
    n = len(Xs)
    Zs = [string[1] for string in Zs_operators]
    
    # Generate all combinations of Xs and their corresponding hats
    for r in range(1, n + 1):  # Start from 1 to include single elements
        for combo in combinations(Xs, r):
            hat = reduce(operator.xor, combo)
            all_x_operators.append([combo, hat])
    
    all_costs = {}

    for used_Xs, X_combos in all_x_operators:
        total_cost = 0
        for Z in Zs:
            cost = ncnot(X_combos | Z)
            total_cost += cost
        
        all_costs[used_Xs] = total_cost
    
    # Find the best combination of Xs that minimizes the cost
    best_Xs = []
    best_cost = 0
    covered = set()
    required = set(Xs)
    maybe_later = []
    
    while len(best_Xs) < n:
        lowest_cost = min(all_costs.values())
        keys = [k for k, v in all_costs.items() if v == lowest_cost]
        
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

    best_Xs_reduced = [reduce(operator.xor, x) for x in best_Xs]

    return best_Xs_reduced, best_cost

if __name__ == '__main__':
    results = find_best_cost([0b0010, 0b0110, 0b1000, 0b1001, 0b1111, 0b0000], [(1, 0b0010), (1, 0b0110), (1, 0b1000), (1, 0b1010), (1, 0b1100), (1, 0b1110)])

    print("Best combo of Xs (heuristic):", results[0],"\nBest cost (heuristic):", results[2])#, "\nBest combo of Xs (exact):", results[2], "\nBest cost (exact):", results[3])
    print("Best combo of Xs reduced (heuristic):", results[1])
    
    # family_of_valid_graphs_1 = {0b0010 : [(0, 1), (2, 7), (3, 9), (4, 13), (5, 10), (6, 8), (11, 14)],
    # 0b0111 : [(0, 2), (1, 7), (3, 4), (5, 14), (6, 12), (9, 13), (10, 11)],
    # 0b1010 : [(0, 3), (1, 9), (2, 4), (6, 11), (7, 13), (8, 14), (10, 12)],
    # 0b1101 : [(0, 4), (1, 13), (2, 3), (5, 8), (6, 10), (7, 9), (11, 12)],
    # 0b1110 : [(0, 5), (1, 10), (2, 14), (4, 8), (6, 13), (7, 11), (9, 12)],
    # 0b0001 : [(0, 6), (1, 8), (2, 12), (3, 11), (4, 10), (5, 13), (9, 14)],
    # 0b0101 : [(0, 7), (1, 2), (3, 13), (4, 9), (5, 11), (8, 12), (10, 14)],
    # 0b0011 : [(0, 8), (1, 6), (3, 14), (4, 5), (7, 12), (9, 11), (10, 13)],
    # 0b1000 : [(0, 9), (1, 3), (2, 13), (4, 7), (5, 12), (6, 14), (8, 11)],
    # 0b1100 : [(0, 10), (1, 5), (2, 11), (3, 12), (4, 6), (7, 14), (8, 13)],
    # 0b1011 : [(0, 11), (1, 14), (2, 10), (3, 6), (4, 12), (5, 7), (8, 9)],
    # 0b0110 : [(0, 12), (2, 6), (3, 10), (4, 11), (5, 9), (7, 8), (13, 14)],
    # 0b1111 : [(0, 13), (1, 4), (2, 9), (3, 7), (5, 6), (8, 10), (12, 14)],
    # 0b1001 : [(0, 14), (1, 11), (2, 5), (3, 8), (6, 9), (7, 10), (12, 13)],
    # 0b0100 : [(1, 12), (2, 8), (3, 5), (4, 14), (6, 7), (9, 10), (11, 13)]}
    # operators_1 = [0b1001, 0b0110, 0b111]
    # nodes_1 = (0, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 14)

    # operators_2 = [0b0010, 0b0110, 0b1000, 0b1010, 0b1100, 0b1110]
    # nodes_2 = (2, 6, 7, 8)
    
    # suborbits = split_into_suborbits(family_of_valid_graphs=family_of_valid_graphs_1, operators=operators_2, nodes=nodes_2)
    # #print("suborbits: ", suborbits)
