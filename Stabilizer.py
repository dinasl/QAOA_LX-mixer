import numpy as np
import itertools
import utils
from functools import reduce
import operator
#TODO update edges
class Stabilizer:
    def __init__(self, familiy_of_graphs, B, n):
        """
        Note: functions need to be called in specific order. get_orbits (Dina's method), then compute_minimal_generating_set
        
        Attributes:
            family_of_graphs (...):
            B ([[]]): Feasible set of bitstrings (int representations) from the computational basis that is an orbit.
            n (int): number of qubits
            orbits ([[ ]]):
            edges ([[[,]]]):
            minimal_generating_set ([[()]]):
            projector ([[[,]]]):

        """
        self.family_of_graphs = familiy_of_graphs
        self.B = B #TODO misledende å kalle de B? er jo egt subsets av B
        self.n = n
        self.orbits = []
        self.nodes_stabilized = []
        self.minimal_generating_sets = []
        self.projectors = []
    
    #helping function to test the code before adding the find_orbit_function
    def set_orbits(self,orbits):
        self.orbits = orbits

    def check_if_orbit(self):
        """
        Checks if the set B is an orbit of a stabilizer group.

        Args:
            B (list[int]): Feasible set of bitstrings (int representations) from the computational basis. 
            
        Returns:
            bool: True if B is an orbit of a stabilizer group, False otherwise.
        """
        #This loop finds if there is an orbit when B=2^n states, iterates over all the possible X's and states, brute force
        #not going to be necessary method
        for B in self.B:
            B_set = set(B)
            tried_X = set()
            for i in range(len(B)-1):
                X_1 = B[i]^B[i+1] #making X_1 = z1 ^ z2
                if X_1 in tried_X:
                    continue
                tried_X.add(X_1)
                    
                #checks if it maps all states to another
                if not all((X_1 ^ state) in B_set for state in B):
                    self.orbits.append([None])
                    break
            self.orbits.append(list(tried_X))
       
    def compute_minimal_generating_sets(self):
        """
        Computes the minimal generating set of a stabilizer group that contains the orbit B.
        
        Args:
            orbits (list[int]): List of x-strings that is the orbit as binary strings (int represenation)
        
        Returns:
            list[int]: List of Pauli strings (int representation) that form the minimal generating set.
        """
        for index, orbit in enumerate(self.orbits):
            #use seed B[0] to get G0 which is on the form G0 = {(+-1, ZII...), ...} where the z-string is on binary (int) form and Z is represented by 1 and I by 0
            #found the seed B[0] from the 0th element of the orbit we are looking at.
            B[0] = self.B[index][0]
            G0 = [((-1 if (B[0] >> (self.n - 1 - i)) & 1 else 1), 1 << (self.n - 1 - i)) for i in range(self.n)]
            
            #iteration process for algoritm 1
            for x_string in orbit:
                G0_elements = [t[1] for t in G0]    #selects all of the elements of G that is a z-string (without +-1)
                G0_signs = [t[0] for t in G0]       #selects the +-1 value

                #is a string that checks if X and Z work on the same qubit for a x-string with all z-strings. Ex: 0100 means X and Z both work on qubit 2 
                commutation_string = [x_string & z_string for z_string in G0_elements]
                
                I_c = []
                I_d = []
                for index, j in enumerate(commutation_string):      #iterates over the elements (binary strings)
                    parity_of_string = utils.parity(j)                    #checks the parity of each string
                    if parity_of_string == 0:
                        I_c.append(index) #appends the position of the commuting string
                    else:
                        I_d.append(index) #appends the position of the anti-commuting string
                #the number of elements that needs to be 
                elements_included = len(G0_elements) - len(I_c) - 1
                
                
                I_d_2 = list(itertools.islice(itertools.combinations(I_d, 2), elements_included))
                I_d_2_Z = [(G0_signs[I_d_2[i][0]]*G0_signs[I_d_2[i][1]],G0_elements[I_d_2[i][0]]|G0_elements[I_d_2[i][1]]) for i in range(elements_included)]
                
                #creates a list of tuples (+-1, Z-string) for commuting pairs  
                I_c_Z = [(G0_signs[i], G0_elements[i]) for i in I_c]

                G_new = I_c_Z + I_d_2_Z
                G0 = G_new
                print(G0)
            
            #finds the final minimal generating set and adds it to the list of minimal generating sets
            final_minimal_generating_set_1_orbit = list(G0)
            self.minimal_generating_sets.append(final_minimal_generating_set_1_orbit)

    def compute_projector_stabilizers(self, restricted = False):
        """
        Computes the restricted projectors using the stabilizer formalism approach.
        
        Args:
            minimal generating set of stabilizer_group (list(tuples[int])): List of Pauli strings (int representation) that form the stabilizer group.
            B (list[int]): Feasible set of bitstrings (int representations) from the computational basis.  
        
        Returns:
            ??? : The restricted projector in the form of a (???, vector) or other suitable representation.
        """
        #TODO remember to normalize the projector
        if restricted:
            for j in range(self.minimal_generating_sets):
                #dimension and number of rows and columns of the matrix
                B_not_stabilized = [x for x in self.B[j] if x not in self.nodes_stabilized[j]]
                minimal_generating_set = self.minimal_generating_sets[j]
                
                matrix = np.zeros((len(B_not_stabilized), len(minimal_generating_set)))
                
                #finding elements for matrix
                for i in B_not_stabilized:
                    matrix_row_binary = [i & stabilizer for stabilizer in minimal_generating_set[:][0]]     #
                    matrix_row = []                                                                         #TODO update this to numpy array so that it is faster
                    for index, string in enumerate(matrix_row_binary):
                        parity_string = utils.parity(string)
                        sign = parity_string*minimal_generating_set[index][0]
                        
                        matrix_row.append(sign)

                    matrix[i] = matrix_row
                print(matrix)
        else: #ignoring global phase?
            #TODO it is supposed to be absolutely all stabilizers with all combinations and the identity?
            projectors = []
            for minimal_generating_set in self.minimal_generating_sets:
                #all possible combinations
                #minimal_generating_set = [(1, 13),(1, 7),(-1, 11)] #tried different minimal_generating_set
                k = len(minimal_generating_set)
                projector = []

                signs, z_strings = zip(*minimal_generating_set)
                signs = np.array(signs)
                z_strings = np.array(z_strings)
                
                # Get all binary combinations (2^k × k)
                all_choices = np.array(list(itertools.product([0, 1], repeat=k)))  # shape (2^k, k)

                for choice in all_choices:
                    # Combine signs of selected generators
                    choice = np.asarray(choice, dtype=bool)
                    selected_signs = signs[choice]
                    total_sign = np.prod(selected_signs) if len(selected_signs) > 0 else 1
                    
                    # Combine Pauli strings using XOR
                    selected_zs = z_strings[choice]
                    if len(selected_zs) == 0:
                        combined_z = 0  # identity
                    else:
                        #reduce applies the function ^ iteratively, and 0 is the inital value
                        combined_z = reduce(lambda a, b: a ^ b, selected_zs, 0)

                    projector.append((total_sign, combined_z))
                
                projectors.append(list(projector))
        
        self.projectors = projectors
B = [0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101]
#compute_minimal_generating_set(B, 4)
B1 = [0b11101, 0b01010, 0b10011, 0b00110]
G = [(-1, 0b00010), (-1, 0b00001), (-1, 0b11000), (1, 0b01100)]
#compute_restricted_projector_stabilizer(B1, 5)

stabilizer = Stabilizer(familiy_of_graphs=None, B=[B], n=4)
stabilizer.check_if_orbit()
stabilizer.compute_minimal_generating_sets()
stabilizer.compute_projector_stabilizers()






    
#def compute_restricted_projector_linalg(stabilizer_group, B):
"""
    Computes the restricted projector using linear algebra approach.
    
    Args:
        stabilizer_group (list[int]): List of Pauli strings (int representation) that form the stabilizer group.
        B (list[int]): Feasible set of bitstrings (int representations) from the computational basis.
    
    Returns:
        ??? : The restricted projector in the form of a (???, vector) or other suitable representation.
"""