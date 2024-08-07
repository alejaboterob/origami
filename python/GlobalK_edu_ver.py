import numpy as np
from scipy.sparse import csr_matrix

from BarKe import BarKe
from FoldKe import FoldKe

def GlobalK_edu_ver(Ui, Node, truss, angles):
    '''The function `GlobalK_edu_ver` calculates the global stiffness matrix and internal forces for a
    truss structure with bars, bends, and folds.
    
    Parameters
    ----------
    Ui
        It seems like the description of the parameters is incomplete. Could you please provide more
    information about the parameters Ui, Node, truss, and angles so that I can better understand the
    function `GlobalK_edu_ver` and assist you further?
    Node
        The `GlobalK_edu_ver` function seems to be calculating the global stiffness matrix and internal
    forces for a truss structure based on the input parameters. Here is a brief overview of the
    parameters used in the function:
    truss
        It looks like the code you provided is a function named `GlobalK_edu_ver` that seems to be related
    to structural analysis using finite element methods. The function calculates the global stiffness
    matrix `K` and the internal forces `IF` for a truss structure based on the input parameters `Ui
    angles
        The `angles` parameter seems to contain information related to bending and folding operations in a
    truss structure. It includes keys such as 'bend', 'fold', 'kpb', 'pb0', 'kpf', 'pf0', and 'CM'.
    These keys likely represent different aspects of
    
    Returns
    -------
        The function `GlobalK_edu_ver` returns two variables `IF` and `K`.
    
    '''
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
        eDof = np.array([np.arange(0, 3) + (truss['Bars'][bel, 0])*3, np.arange(0, 3) + (truss['Bars'][bel, 1])*3], dtype="int").ravel()
        _, Rbe, Kbe = BarKe(Ui[eDof], csr_matrix(truss['B'])[bel, eDof], truss['L'][bel], truss['CM'], truss['A'][bel])
        IFb[eDof, :] = (IFb[eDof, :].T + Rbe).T
        I = np.repeat(eDof, 6).reshape(6, 6).T
        J = I.copy().T
        indi[36*bel:36*(bel+1)] = I.ravel()
        indj[36*bel:36*(bel+1)] = J.ravel()
        kentry[36*bel:36*(bel+1)] = Kbe.ravel()

    Kb = csr_matrix((kentry, (indi, indj)), shape=(3*Nn, 3*Nn), dtype=np.float64)

    # from scipy.io import savemat
    # variables = {
    #     'Kb2': Kb,	
    #     'IFb2': IFb,
    #     'kentry2': kentry,
    #     'indi2': indi,	
    #     'indj2': indj,
    #     'Rbe2': Rbe,
    #     'Kbe2': Kbe,
    # }
    # # Save the variables to a .mat file
    # savemat('variables.mat', variables)

    indi = np.zeros(144*angles['bend'].shape[0])
    indj = indi.copy()
    kentry = indi.copy()
    Lbend = truss['L'][:angles['bend'].shape[0]]
    for d_el in range(angles['bend'].shape[0]):
        eDof = np.array([3*angles['bend'][d_el, :], 3*angles['bend'][d_el, :]+1, 3*angles['bend'][d_el, :]+2]).T.ravel()
        bend = angles['bend'][d_el, :]
        _, Rpe, Kpe = FoldKe(Nodenw, bend, angles['kpb'][d_el], angles['pb0'][d_el], Lbend[d_el], angles['CM'])
        IFp[eDof, :] = (IFp[eDof, :].T + Rpe).T
        I = np.repeat(eDof, 12).reshape(12, 12).T
        J = I.copy().T
        indi[144*d_el:144*(d_el+1)] = I.ravel()
        indj[144*d_el:144*(d_el+1)] = J.ravel()
        kentry[144*d_el:144*(d_el+1)] = Kpe.ravel()

    Kbd = csr_matrix((kentry, (indi, indj)), shape=(3*Nn, 3*Nn), dtype=np.float64)  

    indi = np.zeros(144*angles['fold'].shape[0])
    indj = indi.copy()
    kentry = indi.copy()
    Lfold = truss['L'][angles['bend'].shape[0]]
    for fel in range(angles['fold'].shape[0]):
        eDof = np.array([3*angles['fold'][fel, :], 3*angles['fold'][fel, :]+1, 3*angles['fold'][fel, :]+2]).T.ravel()
        fold = angles['fold'][fel, :]
        _, Rpe, Kpe = FoldKe(Nodenw, fold, angles['kpf'], angles['pf0'][fel], Lfold, angles['CM'])   # Lfold[fel]
        
        IFp[eDof, :] = (IFp[eDof, :].T + Rpe).T

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


