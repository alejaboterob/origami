import numpy as np
from scipy.sparse import coo_matrix

def dirc3d(Node, Ele):
    Ne = Ele.shape[0]
    Nn = Node.shape[0]
    D = np.column_stack([
        Node[Ele[:, 1], 0] - Node[Ele[:, 0], 0], 
        Node[Ele[:, 1], 1] - Node[Ele[:, 0], 1], 
        Node[Ele[:, 1], 2] - Node[Ele[:, 0], 2]
    ])
    L = np.sqrt(D[:, 0]**2 + D[:, 1]**2 + D[:, 2]**2)
    D = np.column_stack([D[:, 0] / L, D[:, 1] / L, D[:, 2] / L])
    # Construct the sparse matrix B
    row = np.repeat(np.arange(Ne), 6).flatten('F')  
    col = np.vstack([3*(Ele[:, 0]+1)-3, 3*(Ele[:, 0]+1)-2, 3*(Ele[:, 0]+1)-1,
                             3*(Ele[:, 1]+1)-3, 3*(Ele[:, 1]+1)-2, 3*(Ele[:, 1]+1)-1]).flatten('F')  # Column-major flattening for MATLAB-like behavior
    # data = np.hstack([D, -D]).T.flatten('F') # Data for matrix entries
    # B = coo_matrix((data, (row, col)), shape=(Ne, 3 * Nn))    
    # data = np.concatenate((D.ravel(), -D.ravel()))
    data = np.hstack([D, -D]).T.flatten('F')
    B = coo_matrix((data, (row, col)), shape=(Ne, 3 * Nn))
    B = -B
    return B, L

# import numpy as np
# from scipy.sparse import csr_matrix

# def dirc3d(Node, Ele):
#     Ne = Ele.shape[0]
#     Nn = Node.shape[0]
#     D = np.column_stack([
#         Node[Ele[:, 1], 0] - Node[Ele[:, 0], 0], 
#         Node[Ele[:, 1], 1] - Node[Ele[:, 0], 1], 
#         Node[Ele[:, 1], 2] - Node[Ele[:, 0], 2]
#     ])
#     L = np.sqrt(D[:, 0]**2 + D[:, 1]**2 + D[:, 2]**2)
#     D = np.column_stack([D[:, 0] / L, D[:, 1] / L, D[:, 2] / L])
#     row_indices = np.repeat(np.arange(Ne), 6)
#     col_indices = np.hstack([
#         3 * Ele[:, 0] - 2, 3 * Ele[:, 0] - 1, 3 * Ele[:, 0], 
#         3 * Ele[:, 1] - 2, 3 * Ele[:, 1] - 1, 3 * Ele[:, 1]
#     ]).ravel()
#     data = np.hstack([D, -D]).ravel()
#     B = csr_matrix((data, (row_indices, col_indices)), shape=(Ne, 3 * Nn))
#     B = -B
#     return B, L