# finite difference coefficients algorithm

# import relevant modules
import numpy as np
import scipy.linalg as sl
from scipy.special import factorial


# given a derivative order M, expansion point x0 and grid points x
def fd_coefficients(M, x0, x):
    # initialise matrix of size M+1 x M+1
    matrix = np.zeros((M+1, M+1))

    # initialising a for loop to iterate through each row
    for i in range(M+1):
        # inner loop to iterate through each column within the row
        for j in range(M+1):
            # update the respective element with the Taylor series coefficient
            matrix[i, j] = ((x[i] - x0)**j) / factorial(j)

    # invert matrix to retrieve weights
    coefficients = sl.inv(matrix)

    return coefficients


# example use
m, x0, x = 4, 0, [-2, -1, 0, 1, 2]
weights = fd_coefficients(m, x0, x)
print(weights)
