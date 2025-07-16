import numpy as np
import itertools
import utils
from functools import reduce
from sympy import Matrix, GF

class Stabilizer:
    def __init__(self, B, n, orbit_dictionary):
        """
        Note: functions need to be called in specific order. get_orbits (Dina's method), then compute_minimal_generating_set
        
        Attributes:
            family_of_graphs (...):
            B ([[]]): Feasible set of bitstrings (int representations) from the computational basis that is an orbit.
            n (int): number of qubits
            orbits ([[ ]]):
            valid_states ([[[,]]]): The states that is defined in the original problem as valid. Does not need to be an orbit
            minimal_generating_set ([[()]]):
            projector ([[[,]]]):

        """
        self.B = B 
        self.n = n
        self.orbit_dictionary = orbit_dictionary

    
    #helping function to test the code before adding the find_orbit_function
    def print_values(self):
        return "orbits: " + str(self.orbits) + "\nMin gen set: " + str(self.minimal_generating_sets) + "\nProjectors: " + str(self.projectors)

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
            orbits_new (Dict{sort_tuple: Class}): Dictionary with nodes as keys (represented in tuples) and the class Orbit as values (containing lists of values)
            B (List[int]): 
        
        Returns:
            list[int]: List of Pauli strings (int representation) that form the minimal generating set.
        """
        for nodes, orbit in self.orbit_dictionary.items():
            #use seed to get G0 which is on the form G0 = {(+-1, ZII...), ...} where the z-string is on binary (int) form and Z is represented by 1 and I by 0
            #found the seed from the 0th element of the tuple, which corresponds to a state saved as a bin int in self.B 
            #TODO reducing the orbit like this does NOT work (makes all the projectors equal????), 
            
            seed = self.B[nodes[0]]
            G0 = [((-1 if (seed >> (self.n - 1 - i)) & 1 else 1), 1 << (self.n - 1 - i)) for i in range(self.n)]
            
            #iteration process for algoritm 1
            for x_string in orbit.Xs: #TODO reduced_orbit_x:, orbit.Xs:
                G0_elements = [t[1] for t in G0]    #selects all of the elements of G that is a z-string (without +-1)
                G0_signs = [t[0] for t in G0]       #selects the +-1 value

                #is a string that checks if X and Z work on the same qubit for a x-string with all z-strings. Ex: 0100 means X and Z both work on qubit 2 
                commutation_string = [x_string & z_string for z_string in G0_elements]

                I_c = []
                I_d = []
                for index, j in enumerate(commutation_string):      #iterates over the elements (binary strings)
                    parity_of_string = utils.parity(j)                    #checks the parity of each string
                    if parity_of_string == 1:
                        I_c.append(index) #appends the position of the commuting string
                    else:
                        I_d.append(index) #appends the position of the anti-commuting string

                #the number of elements that needs to be included
                if len(I_d) > 1: #To be able to combine anti-commuting stabilizer, there needs to be at least 2 anti-commuting stabilizers
                    elements_included = len(G0_elements) - len(I_c) - 1
                
                    I_d_2 = list(itertools.islice(itertools.combinations(I_d, 2), elements_included))
                    I_d_2_Z = [(G0_signs[I_d_2[i][0]]*G0_signs[I_d_2[i][1]],G0_elements[I_d_2[i][0]]|G0_elements[I_d_2[i][1]]) for i in range(elements_included)]
                else:
                    #TODO what happens if there is only 1 anti-commuting stabilizer? -> just add it to the list of commuting stabilizers?
                    I_d_2_Z = [] 
                
                #creates a list of tuples (+-1, Z-string) for commuting pairs  
                I_c_Z = [(G0_signs[i], G0_elements[i]) for i in I_c]


                G_new = I_c_Z + I_d_2_Z
                G0 = G_new
            
            #finds the final minimal generating set and adds it to the list of minimal generating sets
            final_minimal_generating_set_1_orbit = G0 #removed list from list(G0) since it is already a list of tuples
            
            self.orbit_dictionary[nodes].Zs = final_minimal_generating_set_1_orbit
            #print("Minimal generating set for orbit :", nodes, "\nis: ", final_minimal_generating_set_1_orbit)

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
                B_not_stabilized = [x for x in self.valid_states[j] if x not in self.B[j]] #TODO denne blir feil siden B bare er nodene i orbiten, må ha et overordnet set og skrive self.overordnet_set - self.B
                minimal_generating_set = self.minimal_generating_sets[j]
                
                #TODO flip the rows and columns -> better later?
                matrix = np.zeros((len(B_not_stabilized), len(minimal_generating_set)))

                #TODO make directly with 0s and 1s? switch parity back in utils, and do a True/False on minimal_generating_set[index][0]
                for i in B_not_stabilized:
                    #& operator on z-string (stabilizer) and state -> 1s if Z works on 1 and take parity of that
                    #-> -1 if anticommuting and 1 if commuting -> multiply by sign of the stabilizer
                    matrix_row = [utils.parity(i & stabilizer)*minimal_generating_set[index][0] for index, stabilizer in minimal_generating_set[:][1]]     #TODO check that minimal_generating_set[:][1] gives z-string
                    matrix[i] = matrix_row
                
                #convert into matrix of 0s and 1s where 1 -> 0, -1 -> 1
                binary_matrix = (matrix == -1).astype(int)

                #transposes matrix to sove it more easily
                matrix_transposed = Matrix(binary_matrix.T.tolist())

                #creating the vector we want to find (all -1s in the article -> 1s in the matrix)
                target = Matrix([1] * matrix.shape[0])

                #finding the solution
                sol = matrix_transposed.solve_least_squares(target, domain=GF(2))
                
                if sol is not None:
                    indices = [i for i, val in enumerate(sol) if val == 1] # subset of column indices

        else: #TODO ignoring global phase?
            for orbit in self.orbit_dictionary.values():
                
                minimal_generating_set = orbit.Zs
                
                #all possible combinations
                k = len(minimal_generating_set)
                projector = []
                
                # TODO Check if minimal_generating_set is empty, might not be necessary if we ensure an actual orbit (now there can be 12 nodes)
                try: 
                    signs, z_strings = zip(*minimal_generating_set)
                    signs = np.array(signs)
                    z_strings = np.array(z_strings)
                except ValueError:
                    #print("No minimal generating set found for orbit:", orbit)
                    #print("No projectors are needed.")
                    continue
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
                        # Reduce applies the function ^ iteratively, and 0 is the inital value
                        combined_z = reduce(lambda a, b: a ^ b, selected_zs, 0)

                    projector.append((int(total_sign), int(combined_z)))
            
    
                # Updating so that we disregard the minimal generating sets and only keep the projectors 
                orbit.Zs = projector #changed to projector from projectors

                #print("Projector for a given orbit: ", orbit.Zs)

if __name__ == '__main__':
    # TODO does not run, need to make orbit_dictionary first
    B = [0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101]
    orbit_dictionary = {"hei":"hei"}
    #compute_minimal_generating_set(B, 4)
    B1 = [0b11101, 0b01010, 0b10011, 0b00110]
    G = [(-1, 0b00010), (-1, 0b00001), (-1, 0b11000), (1, 0b01100)]
    #compute_restricted_projector_stabilizer(B1, 5)

    stabilizer = Stabilizer(B=B, n=4, orbit_dictionary={})
    #stabilizer.check_if_orbit()
    stabilizer.compute_minimal_generating_sets()
    #stabilizer.compute_projector_stabilizers()

    stabilizer2 = Stabilizer(B=[[0b11000, 0b00100, 0b01101, 0b10001]], n=5)
    stabilizer2.check_if_orbit()
    stabilizer2.compute_minimal_generating_sets()
    stabilizer2.compute_projector_stabilizers()
    print(stabilizer2.print_values())
    



    
#def compute_restricted_projector_linalg(stabilizer_group, B):
"""
    Computes the restricted projector using linear algebra approach.
    
    Args:
        stabilizer_group (list[int]): List of Pauli strings (int representation) that form the stabilizer group.
        B (list[int]): Feasible set of bitstrings (int representations) from the computational basis.
    
    Returns:
        ??? : The restricted projector in the form of a (???, vector) or other suitable representation.
"""