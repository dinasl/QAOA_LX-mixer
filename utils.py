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