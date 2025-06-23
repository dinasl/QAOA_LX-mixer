# Imports

def check_if_orbit(B):
    """
    Checks if the set B is an orbit of a stabilizer group.

    Args:
        B (list[int]): Feasible set of bitstrings (int representations) from the computational basis. 
        
    Returns:
        bool: True if B is an orbit of a stabilizer group, False otherwise.
    """
    #This loop finds if there is an orbit when B=2^n states, iterates over all the possible X's and states
    B_set = set(B)
    tried_X = set()

    for i in range(len(B)-1):
        X_1 = B[i]^B[i+1] #making X_1 = z1 ^ z2
        if X_1 in tried_X:
            continue
        tried_X.add(X_1)
        
        #checks if it maps all states to another
        if not all((X_1 ^ state) in B_set for state in B):
            return False
    
    return True
    #could also make a graph with nodes B, and add edges with X's...(?)
    """
    #brute force:
    for i in range(len(B)-1):
        X_1 = B[i]^B[i+1] #making X_1 = z1 ^ z2
        
        #checks if it maps all, need a way to minimize this one so that it doesn't iterate over the same ones
        for state in B:                     #TODO at the very least we don't need to check the ones that 
            if (X_1 ^ state) not in B:
                return False
    
    return True
    """

def find_orbit_X_string(B):
    pass

def compute_minimal_generating_set(B):
    """
    Computes the minimal generating set of a stabilizer group that contains the orbit B.
    
    Args:
        B (list[int]): Feasible set of bitstrings (int representations) from the computational basis that is an orbit.
    
    Returns:
        list[int]: List of Pauli strings (int representation) that form the minimal generating set.
    """
    pass

#def compute_restricted_projector_linalg(stabilizer_group, B):
    """
    Computes the restricted projector using linear algebra approach.
    
    Args:
        stabilizer_group (list[int]): List of Pauli strings (int representation) that form the stabilizer group.
        B (list[int]): Feasible set of bitstrings (int representations) from the computational basis.
    
    Returns:
        ??? : The restricted projector in the form of a (???, vector) or other suitable representation.
    """

def compute_restricted_projector_stabilizer(minimal_generating_set, B):
    """
    Computes the restricted projector using the stabilizer formalism approach.
    
    Args:
        minimal generating set of stabilizer_group (list[int]): List of Pauli strings (int representation) that form the stabilizer group.
        B (list[int]): Feasible set of bitstrings (int representations) from the computational basis.  
    
    Returns:
        ??? : The restricted projector in the form of a (???, vector) or other suitable representation.
    """
    pass

B = [0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010]
print(check_if_orbit(B))