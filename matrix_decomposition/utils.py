# coding:utf-8
"""
Created by @rh01 on 17-12-15.
author: http://github.com/rh01
"""

from math import sqrt
import numpy as np

def mult_matrix(M, N):
    """Multiply square matrices of same dimension M and N"""

    # Converts N into a list of tuples of columns                                                                                                                                                                                                      
    tuple_N = zip(*N)

    # Nested list comprehension to calculate matrix multiplication                                                                                                                                                                                     
    return [[sum(el_m * el_n for el_m, el_n in zip(row_m, col_n)) for col_n in tuple_N] for row_m in M]

def pivot_matrix(M):
    """Returns the pivoting matrix for M, used in Doolittle's method."""
    m = len(M)

    # Create an identity matrix, with floating point values                                                                                                                                                                                            
    id_mat = [[float(i ==j) for i in xrange(m)] for j in xrange(m)]

    # Rearrange the identity matrix such that the largest element of                                                                                                                                                                                   
    # each column of M is placed on the diagonal of of M                                                                                                                                                                                               
    for j in xrange(m):
        row = max(xrange(j, m), key=lambda i: abs(M[i][j]))
        if j != row:
            # Swap the rows                                                                                                                                                                                                                            
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]

    return id_mat




def trans_matrix(M):
    """Take the transpose of a matrix."""
    n = len(M)
    return [[ M[i][j] for i in range(n)] for j in range(n)]

def norm(x):
    """Return the Euclidean norm of the vector x."""
    return sqrt(sum([x_i**2 for x_i in x]))

def Q_i(Q_min, i, j, k):
    """Construct the Q_t matrix by left-top padding the matrix Q                                                      
    with elements from the identity matrix."""
    if i < k or j < k:
        return float(i == j)
    else:
        return Q_min[i-k][j-k]
    