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