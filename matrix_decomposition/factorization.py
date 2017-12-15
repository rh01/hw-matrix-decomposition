# coding:utf-8
"""
Created by @rh01 on 17-12-15.
author: http://github.com/rh01
"""

from utils import *
   

def LU_decomposition(A, verbose=False):
    """Performs an LU Decomposition of A (which must be square)                                                                                                                                                                                        
    into PA = LU. The function returns P, L and U."""
    n = len(A)

    # Create zero matrices for L and U                                                                                                                                                                                                                 
    L = [[0.0] * n for i in xrange(n)]
    U = [[0.0] * n for i in xrange(n)]

    # Create the pivot matrix P and the multipled matrix PA                                                                                                                                                                                            
    P = pivot_matrix(A)
    PA = mult_matrix(P, A)

    # Perform the LU Decomposition                                                                                                                                                                                                                     
    for j in xrange(n):
        # All diagonal entries of L are set to unity                                                                                                                                                                                                   
        L[j][j] = 1.0

        # LaTeX: u_{ij} = a_{ij} - \sum_{k=1}^{i-1} u_{kj} l_{ik}                                                                                                                                                                                      
        for i in xrange(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in xrange(i))
            U[i][j] = PA[i][j] - s1

        # LaTeX: l_{ij} = \frac{1}{u_{jj}} (a_{ij} - \sum_{k=1}^{j-1} u_{kj} l_{ik} )                                                                                                                                                                  
        for i in xrange(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in xrange(j))
            L[i][j] = (PA[i][j] - s2) / U[j][j]

    
    if verbose :
        print"The LU_decomposition decomposition of A is:"
        print "The L matrix is"
        pprint (L)
        print "The U matrix is"
        pprint (U)

    return L,U


def Gram_Schimidt(A,verbose=False):
    m = len(A)
    # Create zero matrices for L and U                                                                                                                                                                                                                 
    Q = [[0.0] * m for i in xrange(m)]
    R = [[0.0] * m for i in xrange(m)]
    
    R[0][0] = norm(A[:][0])
    #R[0,0]=np.linalg.norm(A[:,0])
    if R[0][0] == 0:
        print"cannot realiaze Gram_Schimidt decomposition"
    else:
        Q[:][0]=[[A[k][0]/R[0][0]] for k in xrange(m)]
    for i in range(1,m):
        Q[:][i]=A[:][i]
        for j in range(i):
            R[j][i]=sum(Q[j][k]*Q[k][i] for k in xrange(m))
            print(R[j][i],Q[:][j])
            Q[:][i]=[[Q[k][i]-R[j][i]*Q[k][j]] for k in xrange(m) ]
        #R[i,i]=np.linalg.norm(Q[:,i])
        R[i][i]=norm(Q[:][i])
        if R[i][i]==0:
            print"cannot realiaze Gram_Schimidt decomposition"
        else:
            Q[:][i] = Q[:][i] / R[i][i]
    if verbose:
        print"The Gram_Schimidt decomposition of A is:"
        print ("The Q matrix is:")
        print (Q)
        print ("The R matrix is:")
        print (R)
    return  Q,R





def householder_reduction(A, verbose=False):
    """Performs a Householder Reflections based QR Decomposition of the                                               
    matrix A. The function returns Q, an orthogonal matrix and R, an                                                  
    upper triangular matrix such that A = QR."""
    n = len(A)

    # Set R equal to A, and create Q as a zero matrix of the same size
    R = A
    Q = [[0.0] * n for i in xrange(n)]

    # The Householder procedure
    for k in range(n-1):  # We don't perform the procedure on a 1x1 matrix, so we reduce the index by 1
        # Create identity matrix of same size as A                                                                    
        I = [[float(i == j) for i in xrange(n)] for j in xrange(n)]

        # Create the vectors x, e and the scalar alpha
        # Python does not have a sgn function, so we use cmp instead
        x = [row[k] for row in R[k:]]
        e = [row[k] for row in I[k:]]
        alpha = -cmp(x[0],0) * norm(x)

        # Using anonymous functions, we create u and v
        u = map(lambda p,q: p + alpha * q, x, e)
        norm_u = norm(u)
        v = map(lambda p: p/norm_u, u)

        # Create the Q minor matrix
        Q_min = [ [float(i==j) - 2.0 * v[i] * v[j] for i in xrange(n-k)] for j in xrange(n-k) ]

        # "Pad out" the Q minor matrix with elements from the identity
        Q_t = [[ Q_i(Q_min,i,j,k) for i in xrange(n)] for j in xrange(n)]

        # If this is the first run through, right multiply by A,
        # else right multiply by Q
        if k == 0:
            Q = Q_t
            R = mult_matrix(Q_t,A)
        else:
            Q = mult_matrix(Q_t,Q)
            R = mult_matrix(Q_t,R)

    # Since Q is defined as the product of transposes of Q_t,
    # we need to take the transpose upon returning it
    

    if verbose:
        print "The householder_reduction of A is:"
        print "The Q matrix is:"
        print Q
        print "The R matrix is:"
        print R
    return trans_matrix(Q), R

def Givens_reduction(A,verbose=False):
    """Performs a Givens Reflections based QR Decomposition of the                                               
    matrix A. The function returns Q, an orthogonal matrix and R, an                                                  
    upper triangular matrix such that A = QR."""
    m = len(A)
    
    # create a identity matrix
    P = [[float(i == j) for i in xrange(m)] for j in xrange(m)]
    R=A
    
    
    for k in range(m):
        for j in range(m-1,k,-1):
            mag = np.sqrt(R[k][k]*R[k][k]+R[j][k]*R[j][k])
            c= R[k][k] / mag
            s= R[j][k] / mag
            I= [[float(i == j) for i in xrange(m)] for j in xrange(m)]
            I[k][k]=c
            I[k][j]=s
            I[j][j]=c
            I[j][k]=-s
            P=mult_matrix(I,P)
            R=mult_matrix(I,R)
    if verbose:
        print"The Givens_reduction of A is:"
        print ("The Q matrix is:")
        print (P)
        print ("The R matrix is:")
        print (R)
    return  trans_matrix(P),R
# A=np.array([[1 , 2 , -3 , 4],[4 , 8 , 12 , -8],[2 , 3 , 2 , 1 ],[-3 , -1 , 1 , -4 ]])
# A = np.array([[1, 19, -34], [-2, -5, 20], [2, 8, 37]])
#A = np.array([[0, -20, -14], [3, 27, -4], [4, 11, -2]])
# A = np.array([[3, 2, 1], [2, -3, 4], [5, 1, -1], [7, 4, 2]])
# A=np.array([[0,0],[0,0]])
if __name__=='__main__':
    #A = np.random.randn(9,9)
    #A = np.random.randint(0,10000,[4,4])
    #A = np.array([[0, -20, -14], [3, 27, -4], [4, 11, -2]])
    A = [[12, -51, 4], [6, 167, -68], [-4, 24, -41]]

    print"the A matrix is:"
    print A
    print "please choose the way to decompose A . "
    print  "Press '1' to realize the PLU_decomposition of A. "
    print  "Press '2' to realize the Gram_Schimidt of A."
    print  "Press '3' to realize the householder_reduction of A."
    print  "Press '4' to realize the Givens_reduction of A."
    #Q = np.empty(A.shape)
    #R = np.empty(A.shape)
    choose_mode=input("please choose the mode:")
    if choose_mode==1:
        L, U=LU_decomposition(A)
        print "The L matrix is"
        print (L)
        print "The U matrix is"
        print (U)

    elif choose_mode==2:
        Q, R  = Gram_Schimidt(A)
        print  np.sum(np.linalg.norm(A - np.dot(Q, R)))

    elif choose_mode==3:
        Q, R = householder_reduction(A)
        print  np.sum(np.linalg.norm(A - np.dot(Q, R)))
    elif choose_mode==4:
        Q,R = Givens_reduction(A)
        print  np.sum(np.linalg.norm(A - np.dot(Q, R)))

    if choose_mode is not 1:
        print ("The Q matrix is:")
        print (Q)
        print ("The R matrix is:")
        print (R)

        print "QR decomposition in numpy is \n", np.linalg.qr(A)
