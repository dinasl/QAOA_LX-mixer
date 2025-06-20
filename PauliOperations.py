from sympy.physics.paulialgebra import Pauli
from sympy.physics.quantum import TensorProduct

from BitString import *
from PauliString import *

def Xoperator(s1, s2):
    """ Take two Bitstrings and return X, such that X*s1=s2
    @param s1 Bitstring
    @param s2 Bitstring
    @return PauliString
    """
    if not (isinstance(s1, BitString) and isinstance(s2, BitString)):
        raise Exception("Input must be of type BitString")
    if s1.scalar!=s2.scalar:
        raise Exception("scalars must be equal")
    X = ""
    for ind in range(len(s1.state)):
        if s1.state[ind] == s2.state[ind]:
            X += "I"
        else:
            X += "X"
    return PauliString(1, X)

# def outer(bitstring):
#     """ Take a Bitstring x=bitstring and calculate |x><x| in the Pauli basis
#     @param bitstring Bitstring
#     @return TensorProduct of Paulis
#     """
#     if not isinstance(bitstring, BitString):
#         raise Exception("Input must be of type BitString")
#     for i in range(len(bitstring)):
#         if bitstring[i]=="0":
#             tmp=1/2*(1+Pauli(3))
#         else:
#             tmp=1/2*(1-Pauli(3))
#         if i == 0:
#             pauli_str=tmp
#         else:
#             pauli_str=TensorProduct(pauli_str,tmp)
#     return pauli_str

# def costPS(PS, H=None):
#     """ Cost of a list of Pauli strings, if H!=None, cost of PS*H
#     @param PS list of PauliString
#     @param H PauliString
#     @return cost=nCNOTs
#     """
#     if not isinstance(PS, list):
#         raise Exception("Input must be a list")

#     ncnot = 0
#     for P in PS:
#         if H !=None:
#             P=H*P
#         ncnot += P.ncnot()
#     return ncnot

