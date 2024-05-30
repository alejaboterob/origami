import numpy as np
from scipy.sparse import csr_matrix

from BarKe import BarKe
from FoldKe import FoldKe

def GlobalK_edu_ver(Ui, Node, truss, angles):
    Nn = Node.shape[0]
    IFb = np.zeros((3*Nn, 1))
    IFp = IFb.copy()
    indi = np.zeros(36*truss['Bars'].shape[0])
    indj = indi.copy()
    kentry = indi.copy()
    Nodenw = Node.copy()
    Nodenw[:, 0] += Ui[::3]
    Nodenw[:, 1] += Ui[1::3]
    Nodenw[:, 2] += Ui[2::3]

    for bel in range(truss['Bars'].shape[0]):
        eDof = np.array([np.arange(0, 3) + (truss['Bars'][bel, 0])*3, np.arange(0, 3) + (truss['Bars'][bel, 1])*3]).ravel()
        _, Rbe, Kbe = BarKe(Ui[eDof], csr_matrix(truss['B'])[bel, eDof], truss['L'][bel], truss['CM'], truss['A'][bel])
        IFb[eDof, :] = (IFb[eDof, :].T + Rbe).T
        I = np.repeat(eDof, 6).reshape(6, 6).T
        J = I.copy().T
        indi[36*bel:36*(bel+1)] = I.ravel()
        indj[36*bel:36*(bel+1)] = J.ravel()
        kentry[36*bel:36*(bel+1)] = Kbe.ravel()

    Kb = csr_matrix((kentry, (indi, indj)), shape=(3*Nn, 3*Nn), dtype=np.float64)

    indi = np.zeros(144*angles['bend'].shape[0])
    indj = indi.copy()
    kentry = indi.copy()
    Lbend = truss['L'][:angles['bend'].shape[0]]
    for d_el in range(angles['bend'].shape[0]):
        eDof = np.array([3*angles['bend'][d_el, :], 3*angles['bend'][d_el, :]+1, 3*angles['bend'][d_el, :]+2]).T.ravel()
        bend = angles['bend'][d_el, :]
        _, Rpe, Kpe = FoldKe(Nodenw, bend, angles['kpb'], angles['pb0'][d_el], Lbend[d_el], angles['CM'])
        IFb[eDof, :] = (IFb[eDof, :].T + Rpe).T
        I = np.repeat(eDof, 12).reshape(12, 12).T
        J = I.copy().T
        indi[144*d_el:144*(d_el+1)] = I.ravel()
        indj[144*d_el:144*(d_el+1)] = J.ravel()
        kentry[144*d_el:144*(d_el+1)] = Kpe.ravel()

    Kbd = csr_matrix((kentry, (indi, indj)), shape=(3*Nn, 3*Nn), dtype=np.float64)  
    if Kbd.size == 0:
        Kbd = np.zeros((3*Nn, 3*Nn))

    indi = np.zeros(144*angles['fold'].shape[0])
    indj = indi.copy()
    kentry = indi.copy()
    Lfold = truss['L'][angles['bend'].shape[0]:]
    for fel in range(angles['fold'].shape[0]):
        eDof = np.array([3*angles['fold'][fel, :], 3*angles['fold'][fel, :]+1, 3*angles['fold'][fel, :]+2]).T.ravel()
        fold = angles['fold'][fel, :]
        _, Rpe, Kpe = FoldKe(Nodenw, fold, angles['kpf'], angles['pf0'][fel], Lfold[fel], angles['CM'])
        
        IFb[eDof, :] = (IFb[eDof, :].T + Rpe).T

        I = np.repeat(eDof, 12).reshape(12, 12).T
        J = I.copy().T
        indi[144*fel:144*(fel+1)] = I.ravel()
        indj[144*fel:144*(fel+1)] = J.ravel()
        kentry[144*fel:144*(fel+1)] = Kpe.ravel()

    Kfd = csr_matrix((kentry, (indi, indj)), shape=(3*Nn, 3*Nn))

    IF = IFb + IFp
    K = Kb + Kbd + Kfd
    K = (K + K.T) / 2

# from scipy.io import savemat
# variables = {
#     'Kbd2': Kbd,
#     'kentry2': kentry,
#     'indi2': indi,
#     'indj2': indj,

# }
# # Save the variables to a .mat file
# savemat('Kbd.mat', variables)


# from scipy.io import savemat
# variables = {
#     'Kb2': Kb,
#     'Kbd2': Kbd,
#     'Kfd2': Kfd,
#     'IFb2': IFb,
#     'IFp2': IFp,
#     'IF2': IF,
#     'K2': K,
# }
# # Save the variables to a .mat file
# savemat('variables.mat', variables)

    return IF, K


