from sympy import *
from sympy.physics.paulialgebra import Pauli
from sympy.physics.quantum import TensorProduct

from BitString import *

class PauliString:

    def __init__(self, scalar, P):
        if not (isinstance(scalar, complex) or isinstance(scalar, float) or isinstance(scalar, int)):
            raise Exception("Scalar must be type complex, float or int, but received ", type(scalar))
        if isinstance(P, TensorProduct) or isinstance(P, Pauli):
            self.P = PauliString.__TensorProductPaulis2String(P)
        elif isinstance(P, str):
            if not all(c in 'IXYZ' for c in P):
                raise Exception("Pauli string must contain only I, X, Y, Z. Actual: ", P)
            self.P = P
        else:
            raise Exception("State must be Pauli, TensorProduct or str, but received ", type(P))
        self.scalar = scalar

    def __eq__(self, other):
        if isinstance(other, PauliString):
            return self.scalar.real == other.scalar.real and self.scalar.imag == other.scalar.imag and self.P == other.P
        return False

    def __str__(self):
        return str(self.scalar) + "*" + self.P

    def __repr__(self):
        return repr(self.scalar) + "*" + self.P

    def __hash__(self):
        return hash(str(self))

    @staticmethod
    def __PaulitoString(item):
        if item==1:
            return "I"
        elif item==Pauli(3):
            return "Z"
        elif item==Pauli(1):
            return "X"
        elif item==Pauli(2):
            return "Y"

    @staticmethod
    def __TensorProductPaulis2String(P):
        args=P.args
        tmp=''
        if len(args)==0:
            tmp+=PauliString.__PaulitoString(P)
        else:
            if isinstance(args[0],TensorProduct):
                tmp+=PauliString.__TensorProductPaulis2String(args[0])
            else:
                tmp+=PauliString.__PaulitoString(args[0])
            if isinstance(args[1],TensorProduct):
                tmp+=PauliString.__TensorProductPaulis2String(args[1])
            else:
                tmp+=PauliString.__PaulitoString(args[1])
        return tmp


    def checkreal(self):
        if self.scalar.imag==0:
            self.scalar = self.scalar.real
            return True
        else:
            return False

    def ncnot(self):
        ncnot = len(self.P.replace("I", ""))
        if ncnot > 1:
            return 2 * (ncnot - 1)
        else:
            return 0

    @staticmethod
    def __mulP2P(P1, P2):
        if P1==P2:
            return  1, 'I'
        elif P1=='I':
            return 1, P2
        elif P2=='I':
            return 1, P1
        elif P1=='X':
            if P2=='Y':
                return 1j, 'Z'
            else:## only case left is P2==Z
                return -1j, 'Y'
        elif P1=='Y':
            if P2=='X':
                return -1j, 'Z'
            else:## only case left is P2==Z
                return 1j, 'X'
        elif P1=='Z':
            if P2=='X':
                return 1j, 'Y'
            else:## only case left is P2==Y
                return -1j, 'X'

    @staticmethod
    def __mulP2S(P,s):
        if P=='I':
            return 1, s
        elif P=='X':
            if s=='0':
                return 1, '1'
            else:
                return 1, '0'
        elif P=='Y':
            if s=='0':
                return 1j, '1'
            else:
                return -1j, '0'
        elif P=='Z':
            if s=='0':
                return +1, '0'
            else:
                return -1, '1'

    def __rmul__(self, other):
        if isinstance(other, complex) or isinstance(other, float) or isinstance(other, int):
            return PauliString(other*self.scalar, self.P)
        else:
            raise Exception("Multiplication from the right with "+str(type(other))+" not implemented.")
        
    def __mul__(self, other):
        if isinstance(other,PauliString):
            if len(self.P)!=len(other.P):
                raise Exception("Length of Pauli strings must be equal.")
            val=self.scalar*other.scalar
            PS=''
            for i in range(len(self.P)):   
                f, P = self.__mulP2P(self.P[i], other.P[i])
                val*=f
                PS+=P
            return PauliString(val, PS)
        elif isinstance(other,BitString):
            if len(self.P)!=len(other.state):
                raise Exception("Length of Pauli string and bit string must be equal.")
            val=self.scalar*other.scalar
            bs=''
            for i in range(len(self.P)):   
                f, s = self.__mulP2S(self.P[i], other.state[i])
                val*=f
                bs+=s
            return BitString(val,bs)
        else:
            raise Exception("Multiplication with "+str(type(other))+" not implemented.")

