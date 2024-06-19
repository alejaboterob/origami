import numpy as np
from scipy.sparse import csr_matrix
from findbend import findbend
from findfdbd import findfdbd
from dirc3d import dirc3d
from FoldKe import FoldKe

def PrepareData(Node, Panel, Supp, Load, BarCM, RotSpring, kpf, kpb, Abar):
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
        'kpf': kpf,
        'kpb': kpb,
        'pf0': pf0,
        'pb0': np.full(Bend.shape[0], np.pi),
        'Panel': Panel
    }

    return truss, angles, F
