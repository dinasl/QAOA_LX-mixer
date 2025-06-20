class BitString:

    def __init__(self, scalar, state):
        if not (isinstance(scalar, complex) or isinstance(scalar, float) or isinstance(scalar, int)):
            raise Exception("scalar must be type float.")
        if not isinstance(state, str):
            raise Exception("state must be type str.")
        if not all(c in '01' for c in state):
            raise Exception("state must contain only 0 and 1.")
        self.scalar = scalar
        self.state = state

    def __eq__(self, other):
        if isinstance(other, BitString):
            return self.scalar.real == other.scalar.real and self.scalar.imag == other.scalar.imag and self.state == other.state
        return False

    def __str__(self):
        return str(self.scalar) + "*" + self.state

    def __repr__(self):
        return repr(self.scalar) + "*" + self.state
