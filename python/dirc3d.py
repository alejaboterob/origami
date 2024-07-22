import numpy as np

from scipy.sparse import coo_matrix

def dirc3d(Node, Ele):
    '''The function `dirc3d` calculates the direction cosines of each element in a 3D mesh and constructs a
    sparse matrix B based on the calculated values.

    Parameters
    ----------
    Node
        The `Node` parameter in the `dirc3d` function likely represents the coordinates of the nodes in a
    3D space. It is a numpy array where each row corresponds to the coordinates of a node in 3D space.
    Ele
        The `Ele` parameter in the `dirc3d` function represents the element connectivity array. It is a
    NumPy array where each row defines an element by specifying the indices of the nodes that form the
    element.

    Returns
    -------
        The function `dirc3d` returns two values: 
    1. The sparse matrix `B` which represents the directional cosines of the elements in a 3D mesh.
    2. The array `L` which contains the lengths of the elements in the mesh.

    '''
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