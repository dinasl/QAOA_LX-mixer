is_power_of_two = lambda x: (x > 0) and (x & (x - 1)) == 0

# Perhaps use string representation rather than binary int representation for Pauli strings (for better reusability with projectors)??

def ncnot(P) :
    """
    Calculate the number of CNOT gates required to implement a Pauli string.

    Args:
        P (int): Pauli string represented as a binary integer, where each 1 represents a qubit that is acted upon by a Pauli operator (X or Z).

    Returns:
        int: Number of CNOT gates required to implement the Pauli string.
    """
    ncnot = len(str(P).replace("0", ""))
    if ncnot > 1:
        return 2 * (ncnot - 1)
    else:
        return 0

# def costPS(PS, h=None):
#     """
#     Calculate the cost of a list of Pauli strings.

#     Args:
#         PS (array-like[int]): Array of pauli strings (binary int representation).
#         h (int, optional): A constant multiplier for the cost. Defaults to None.

#     Returns:
#         int: The cost of the Pauli string.
#     """
#     # Is h an int or float?
#     if h!=None: PS*= h
#     ncnots = ncnot(PS)
#     return sum(ncnots)