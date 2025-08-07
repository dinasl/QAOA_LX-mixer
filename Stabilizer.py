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
            orbit_dictionary (Dict{sort_tuple: Class}): Dictionary with nodes as keys (represented in tuples) and the class Orbit as values (containing lists of values)

        """
        self.B = B 
        self.n = n
        self.orbit_dictionary = orbit_dictionary

    def compute_minimal_generating_sets(self):
        """
        Computes the minimal generating set of a stabilizer group that contains the orbit B and updates the orbit dataclass instance with the minimal generating set for Zs.
        """
        for nodes, orbit in self.orbit_dictionary.items():

            #find the seed from the 0th element of the tuple, which corresponds to a state saved as a bin int in self.B
            seed = self.B[nodes[0]]
            
            #use seed to get G0 which is on the form G0 = {(+-1, ZII...), ...} where the z-string is on binary (int) form and Z is represented by 1 and I by 0
            G0 = [((-1 if (seed >> (self.n - 1 - i)) & 1 else 1), 1 << (self.n - 1 - i)) for i in range(self.n)]

            #iteration process for algoritm 1
            for x_string in orbit.Xs: #TODO reduced_orbit_x:, orbit.Xs:
                G0_elements = [t[1] for t in G0]    #selects all of the elements of G that is a z-string (without +-1)
                G0_signs = [t[0] for t in G0]       #selects the +-1 value

                #is a string that checks if X and Z work on the same qubit for a x-string with all z-strings. Ex: 0100 means X and Z both work on qubit 2 
                commutation_string = [x_string & z_string for z_string in G0_elements]

                #I_c and I_d are lists that will contain the indices of the commuting and anti-commuting stabilizers respectively
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
                    I_d_2_Z = [(G0_signs[I_d_2[i][0]]*G0_signs[I_d_2[i][1]],G0_elements[I_d_2[i][0]]^G0_elements[I_d_2[i][1]]) for i in range(elements_included)]
                else:
                    # If there are 0 or 1 anti-commuting stabilizers, we just take the commuting ones
                    I_d_2_Z = [] 
                
                #creates a list of tuples (+-1, Z-string) for commuting pairs  
                I_c_Z = [(G0_signs[i], G0_elements[i]) for i in I_c]
                
                #the new minimal generating set is the combination of the commuting and anti-commuting stabilizers
                G_new = I_c_Z + I_d_2_Z
                
                #updates G0 to the new minimal generating set
                G0 = G_new
            
            #finds the final minimal generating set and adds it to the list of minimal generating sets
            final_minimal_generating_set_1_orbit = G0 
            
            #updates the orbit dataclass instance with the final minimal generating set
            self.orbit_dictionary[nodes].Zs = final_minimal_generating_set_1_orbit

    def compute_projector_stabilizers(self, restricted = False):
        """
        Computes the restricted projectors using the stabilizer formalism approach and updates the orbit dataclass instance with the projectors for Zs. 
        """
        for orbit in self.orbit_dictionary.values():
                
            minimal_generating_set = orbit.Zs
                
            #all possible combinations
            k = len(minimal_generating_set)
            projector = []

            #if |B| = 2^n we have an edge case where there are no minimal generating sets and we don't need a projector as B is the entire space. We therefore automatically return the identity projector
            if k == 0:
                orbit.Zs = [(1,0)]  
                return
           
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
                    # Reduce applies the function ^ iteratively, and 0 is the inital value
                    combined_z = reduce(lambda a, b: a ^ b, selected_zs, 0)

                projector.append((int(total_sign), int(combined_z)))

            # Updating the orbit dataclass instance so that we disregard the minimal generating sets and only keep the projectors for Zs
            orbit.Zs = projector 