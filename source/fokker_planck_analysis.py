from utils import *

def is_diagonal_mat(A):
    return np.count_nonzero(A - np.diag(np.diagonal(A))) == np.prod(A.shape)

def fp_analyze(X, V):
    print('X shape:', X.shape, 'V shape:', V.shape)
    sample_num = X.shape[0]
    # calc A
    A_shape = (X[0].reshape([-1, 1]) @ X[0].reshape([1, -1])).shape
    print('A_shape:', A_shape)
    A = np.zeros(A_shape)
    for i in range(sample_num):
        A += X[i].reshape([-1, 1]) @ X[i].reshape([1, -1])
    A /= sample_num
    # check if A is diagonal
    print("A is diagonal:", is_diagonal_mat(A))
    
    # calc B
    B_shape = (X[0].reshape([-1, 1]) @ X[0].reshape([1, -1])).shape
    print('B_shape:', B_shape)
    B = np.zeros(B_shape)
    for i in range(sample_num):
        B += X[i].reshape([-1, 1]) @ V[i].reshape([1, -1])
    B /= sample_num
        

    # calc F from A and B
    A_inv = np.linalg.inv(A)
    F = B @ A_inv
