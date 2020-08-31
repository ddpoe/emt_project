from utils import *
from eigen_analysis import *
def is_diagonal_mat(A):
    return np.count_nonzero(A - np.diag(np.diagonal(A))) == np.prod(A.shape)

def fp_analyze(X, V, need_log_transform=True):
    '''
    scvelo uses nature log normalization
    input X, V are in PCA space and all log transformed version
    '''
    print('[DEBUG] In fp analyze')
    print('X shape:', X.shape, 'V shape:', V.shape)
    sample_num = X.shape[0]
    V = np.copy(V) # do not change original V
            
    # calc A
    A_shape = (X[0].reshape([-1, 1]) @ X[0].reshape([1, -1])).shape
    print('A_shape:', A_shape)
    A = np.zeros(A_shape)
    for i in range(sample_num):
        A += X[i].reshape([-1, 1]) @ X[i].reshape([1, -1])
    A /= sample_num
    # check if A is diagonal
    print("A is diagonal:", is_diagonal_mat(A))
    print('A nonzerorate:', np.count_nonzero(A) / np.prod(A.shape))
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
    analyze_jacob_eigen_complex_plane(F, prefix='originalF')


    # D part
    D = -1/2 * (B + B.T)
    Q = -1/2 * (B - B.T)
    U = A_inv
    
    D_eigenvals, D_eigenvectors = calc_eigen(D)
    D_reals = np.array([num.real for num in D_eigenvals])
    D_eigenvals[D_reals < 0] = 0
    D_ = D_eigenvectors.T @ np.diag(D_eigenvals) @ D_eigenvectors
    F_ = - (D_ + Q) @ U
    analyze_jacob_eigen_complex_plane(F_, prefix='adjustedF')

    F_eigenvals, F_eigenvectors = calc_eigen(F_)
    F_reals = np.array([num.real for num in F_eigenvals])
    print('#F eigenval real part > 0:', np.sum(F_reals > 0))
    
