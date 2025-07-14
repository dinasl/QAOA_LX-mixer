import math
import networkx as nx

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
    
    print("suborbit size: ", suborbit_size)
    print("Number of nodes: ", len(nodes))
    
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
    suborbits = split_into_suborbits(family_of_valid_graphs={0b0010 : [(0, 1), (2, 7), (3, 9), (4, 13), (5, 10), (6, 8), (11, 14)],
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
    0b0100 : [(1, 12), (2, 8), (3, 5), (4, 14), (6, 7), (9, 10), (11, 13)]}, operators=[0b1001, 0b0110, 0b111], nodes=(0, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 14))
    print("suborbits: ", suborbits)