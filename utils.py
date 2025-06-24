is_power_of_two = lambda x: (x > 0) and (x & (x - 1)) == 0

# Use string representation rather than binary int representation for Pauli strings (for better reusability with projectors)

def ncnot(P, h = None) :
    ncnot = len(str(P).replace("I", ""))
    if ncnot > 1:
        return 2 * (ncnot - 1)
    else:
        return 0

def costPS(PS, h=None):
    """
    Calculate the cost of a Pauli string.

    Args:
        PS (array-like[int]): Array of pauli strings (binary int representation).
        h (int, optional): A constant multiplier for the cost. Defaults to None.

    Returns:
        int: The cost of the Pauli string.
    """
    # Is h an int or float?
    if h!=None: PS*= h
    ncnots = ncnot(PS)
    return sum(ncnots)