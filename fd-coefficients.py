# finite difference coefficients algorithm

# import relevant modules
import numpy as np
from scipy.linalg import inv
from scipy.special import factorial
import sympy as sp


# given a derivative order M, expansion point x0 and grid points x
def fd_coefficients(M, x0, x):

    # create an array of the powers from 0 to M
    powers = np.arange(M + 1)

    # use broadcasting to calculate the matrix
    # this removes the need for nested for loops
    matrix = ((x[:, None] - x0) ** powers) / factorial(powers)

    # invert the matrix to retrieve the weights
    coefficients = inv(matrix)

    # convert to a SymPy matrix for output
    coefficients = sp.Matrix(coefficients)

    # rationalize the entries
    coefficients = coefficients.applyfunc(lambda x: sp.nsimplify(x, rational=True))

    # replace small floating-point numbers close to zero with symbolic zero
    coefficients = coefficients.applyfunc(lambda x: 0 if abs(x) < 1e-10 else x)

    return coefficients

# show matrix using pprint
#sp.pprint(coefficients)
