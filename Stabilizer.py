import numpy as np
import itertools
import utils

#TODO make class

class Stabilizer:

    def __init__(self, B, n):
        self.B = B
        self.n = n
        self.orbit = None
        self.minimal_generating_set = None

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

        B_set = set(self.B)
        tried_X = set()
        for i in range(len(self.B)-1):
            X_1 = self.B[i]^self.B[i+1] #making X_1 = z1 ^ z2
            if X_1 in tried_X:
                continue
            tried_X.add(X_1)
                
            #checks if it maps all states to another
            if not all((X_1 ^ state) in B_set for state in self.B):
                return False, None
        
        return True, tried_X
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

    def compute_minimal_generating_set(self):
        """
        Computes the minimal generating set of a stabilizer group that contains the orbit B.
        
        Args:
            B (list[int]): Feasible set of bitstrings (int representations) from the computational basis that is an orbit.
            n (int): number of qubits
            orbit (list[int]): List of x-strings that is the orbit as binary strings (int represenation)
        
        Returns:
            list[int]: List of Pauli strings (int representation) that form the minimal generating set.
        """
        """
        #TODO what if it is not true?
        if check_if_orbit(B)[0]:
            orbit = check_if_orbit(B)[1]
        """
        
        #use seed B[0] to get G0 which is on the form G0 = {(+-1, ZII...), ...} where the z-string is on binary (int) form and Z is represented by 1 and I by 0
        G0 = [((-1 if (self.B[0] >> (self.n - 1 - i)) & 1 else 1), 1 << (self.n - 1 - i)) for i in range(self.n)]

        #iteration process for algoritm 1
        for x_string in self.orbit:
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

            #creates the anti-commuting pairs, it is now a list withing a list: ex: [[2, 3], [2, 4], ...]
            I_d_2 = list(itertools.combinations(I_d, 2)) 
            
            #only iterating over necessary pairs and stores the relevant ones, TODO do elements included above in I_d_2
            elements_included = len(G0_elements) - len(I_c) - 1
            I_d_2_shortened_Z = []

            if len(I_d_2) > 0: #checking that we have anti-commuting pairs #TODO I am not sure this one is necessary anymore...
                for i in range(elements_included): 
                    #creates a tuple with (+-1, ZiZj) for anti-commuting pairs 
                    #+-1 is the signs multiplied and ZiZj is the bitstrings combined
                    I_d_2_shortened_Z.append((G0_signs[I_d_2[i][0]]*G0_signs[I_d_2[i][1]],G0_elements[I_d_2[i][0]]|G0_elements[I_d_2[i][1]]))

            #creates a list of tuples (+-1, Z-string) for commuting pairs  
            I_c_Z = [(G0_signs[i], G0_elements[i]) for i in I_c]

            G_new = I_c_Z + I_d_2_shortened_Z
            G0 = G_new
        
        return G0
                
                
    def compute_restricted_projector_stabilizer(self):
        """
        Computes the restricted projector using the stabilizer formalism approach.
        
        Args:
            minimal generating set of stabilizer_group (list(tuples[int])): List of Pauli strings (int representation) that form the stabilizer group.
            B (list[int]): Feasible set of bitstrings (int representations) from the computational basis.  
        
        Returns:
            ??? : The restricted projector in the form of a (???, vector) or other suitable representation.
        """
        matrix = np.zeros((len(B), len(self.minimal_generating_set[:][1])))
        #finding elements for matrix
        for i in B:
            matrix_row_binary = [i & stabilizer for stabilizer in self.minimal_generating_set] 
            matrix_row = []                                                                          #TODO update this to numpy array so that it is faster
            for index, string in enumerate(matrix_row_binary):
                parity_string = utils.parity(string)
                sign = parity_string*self.minimal_generating_set[index][0]
                
                matrix_row.append(sign)

            matrix[i] = matrix_row
        print(matrix)



B = [0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101]
#compute_minimal_generating_set(B, 4)
B1 = [0b11101, 0b01010, 0b10011, 0b00110]
G = [(-1, 0b00010), (-1, 0b00001), (-1, 0b11000), (1, 0b01100)]
#compute_restricted_projector_stabilizer(B1, 5)








    
#def compute_restricted_projector_linalg(stabilizer_group, B):
"""
    Computes the restricted projector using linear algebra approach.
    
    Args:
        stabilizer_group (list[int]): List of Pauli strings (int representation) that form the stabilizer group.
        B (list[int]): Feasible set of bitstrings (int representations) from the computational basis.
    
    Returns:
        ??? : The restricted projector in the form of a (???, vector) or other suitable representation.
"""