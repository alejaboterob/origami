import numpy as np
from scipy.sparse import csr_matrix
from findbend import findbend
from findfdbd import findfdbd
from dirc3d import dirc3d
from FoldKe import FoldKe

def PrepareData(Node, Panel, Supp, Load, BarCM, RotSpring, kpf, kpb, Abar):
    '''The function `PrepareData` prepares data for a truss structure analysis by finding hinges, defining
    bar elements, calculating forces, and organizing data into dictionaries.
    
    Parameters
    ----------
    Node
        The `Node` parameter in the `PrepareData` function likely represents the nodes or vertices of the
    truss structure. These nodes define the spatial coordinates of the points where the bars and panels
    of the truss are connected.
    Panel
        Panel is a numpy array containing information about the panels in the truss structure. Each row of
    the Panel array represents a panel and contains data such as panel coordinates or indices of nodes
    that form the panel.
    Supp
        The `Supp` parameter in the `PrepareData` function represents the support conditions for the truss
    structure. It contains information about the support nodes and their corresponding constraints. The
    support conditions are used to determine the fixed degrees of freedom in the truss analysis.
    Load
        The `Load` parameter in the `PrepareData` function represents the applied loads on the nodes of the
    truss structure. It is a numpy array where each row represents a different node and the columns
    represent the node index, the applied load in the x-direction, y-direction, and z-direction
    respectively
    BarCM
        BarCM is the center of mass for the bars in the truss structure.
    RotSpring
        RotSpring is a parameter that represents the rotational spring stiffness at folding hinges in the
    truss structure.
    kpf
        It seems like the description of the parameter `kpf` got cut off. Could you please provide more
    information or context about what `kpf` represents or how it is used in the function `PrepareData`?
    kpb
        The parameter `kpb` in the `PrepareData` function represents the rotational stiffness of the
    bending hinges in the truss structure. It is used to calculate the initial bending stiffness for the
    bending hinges.
    Abar
        The parameter `Abar` in the `PrepareData` function represents the cross-sectional area of the bars
    in the truss structure. It is used to define the area of each bar element in the truss.
    
    Returns
    -------
        The function `PrepareData` returns three dictionaries: `truss`, `angles`, and `F`.
    
    '''
    # Find bending hinges
    Bend = findbend(Panel, Node)
    # Find folding hinges and boundaries, return final triangulation
    Fold, Bdry, Trigl = findfdbd(Panel, Bend)
    # Define bar elements
    Bars = np.vstack([Bend[:, :2], Fold[:, :2], Bdry])
    B, L = dirc3d(Node, Bars.astype(int))

    if Supp.shape[0] == 0:
        rs = []
    else:
        rs = np.hstack([Supp[:, 0:1]*3, Supp[:, 0:1]*3+1, Supp[:, 0:1]*3+2]).flatten()
        rs = np.vstack([rs, Supp[:, 1:].flatten()])
        rs = rs[:,rs[1,:] != 0][0]

    if np.isscalar(Abar):
        Abar = np.full(Bars.shape[0], Abar)

    pf0 = np.zeros(Fold.shape[0])
    for i in range(Fold.shape[0]):
        pf0[i], _, _ = FoldKe(Node, Fold[i,:], kpf, 0)

    pb0 = np.zeros(Bend.shape[0])
    for i in range(Bend.shape[0]):
        pb0[i], _, _ = FoldKe(Node, Bend[i, :], kpb, 0)

    m = Node.shape[0]
    F = np.zeros(3 * m)
    indp = Load[:, 0].astype(int)
    F[3 * indp ] = Load[:, 1]
    F[3 * indp +1] = Load[:, 2]
    F[3 * indp +2 ] = Load[:, 3]

    truss = {
        'CM': BarCM,
        'Node': Node,
        'Bars': Bars,
        'Trigl': Trigl,
        'B': B,
        'L': L,
        'FixedDofs': np.unique(rs),
        'A': Abar
    }

    angles = {
        'CM': RotSpring,
        'fold': Fold,
        'bend': Bend,
        'kpf': np.full(Fold.shape[0],kpf),
        'kpb': np.full(Bend.shape[0],kpb),
        'pf0': pf0,
        # 'pb0': np.float64(np.pi),
        'pb0': pb0*0+np.pi,   # OJOOOO 
        'Panel': Panel
    }

    return truss, angles, F
