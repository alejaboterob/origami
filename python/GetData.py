import numpy as np
from scipy.sparse import csr_matrix
from BarKe import BarKe
from FoldKe import FoldKe

def GetData(Ui, Node, truss, angles):
    '''The function `GetData` calculates various properties of bars and angles in a truss structure based
    on input parameters.
    
    Parameters
    ----------
    Ui
        The `Ui` parameter seems to represent the displacement vector of the truss structure. It is used to
    calculate the updated nodal positions in the `GetData` function.
    Node
        The `Node` parameter represents the nodal coordinates of the truss structure. It is a numpy array
    where each row corresponds to a node and the columns represent the x, y, and z coordinates of the
    node in 3D space.
    truss
        The `GetData` function takes four parameters: `Ui`, `Node`, `truss`, and `angles`.
    angles
        The `angles` parameter seems to be a dictionary containing keys such as 'fold' and 'bend'. Each key
    corresponds to an array of angles. The 'fold' key seems to represent folding angles, while the
    'bend' key represents bending angles.
    
    Returns
    -------
        The function `GetData` returns five arrays: `Exbar`, `FdAngle`, `BdAngle`, `LFd`, and `LBd`.
    
    '''
    Exbar = np.zeros((truss['Bars'].shape[0], 1))
    FdAngle = np.zeros((angles['fold'].shape[0], 1))
    BdAngle = np.zeros((angles['bend'].shape[0], 1))
    LFd = FdAngle.copy()
    LBd = BdAngle.copy()
    
    Nodenw = Node.copy()
    Nodenw[:, 0] = Node[:, 0] + Ui[::3]
    Nodenw[:, 1] = Node[:, 1] + Ui[1::3]
    Nodenw[:, 2] = Node[:, 2] + Ui[2::3]
    
    for bel in range(truss['Bars'].shape[0]):
        eDof = np.array([np.arange(0, 3) + truss['Bars'][bel, 0] * 3,
                        np.arange(0, 3) + truss['Bars'][bel, 1] * 3], dtype="int")
        Exbar[bel], Rbe, Kbe= BarKe(Ui[eDof.ravel()], csr_matrix(truss['B'])[bel, eDof.ravel()], truss['L'][bel], truss['CM'], truss['A'][bel])
    
    for d_el in range(angles['bend'].shape[0]):
        bend = angles['bend'][d_el, :]
        BdAngle[d_el], Rhe, Khe = FoldKe(Nodenw, bend, angles['kpb'], angles['pb0'])
        # BdAngle[d_el], Rhe, Khe = FoldKe(Nodenw, bend, angles['kpb'], angles['pb0'][d_el])
        LBd[d_el] = np.linalg.norm(Nodenw[bend[1]] - Nodenw[bend[0]])
    
    for fel in range(angles['fold'].shape[0]):
        fold = angles['fold'][fel, :]
        FdAngle[fel],  Rhe, Khe = FoldKe(Nodenw, fold, angles['kpf'], angles['pf0'][fel])
        LFd[fel] = np.linalg.norm(Nodenw[fold[1]] - Nodenw[fold[0]])
    
    return Exbar, FdAngle, BdAngle, LFd, LBd

