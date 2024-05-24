import numpy as np
from scipy.sparse import csr_matrix

def findfdbd(Panel, bend):
    Nn = max(np.max(panel) for panel in Panel)+1 # sin +1 igual que matlab
    
    # triangularization
    Panelsize = [len(panel) for panel in Panel]
    Ptri = [Panel[i] for i in range(len(Panelsize)) if Panelsize[i] == 3]
    # Triglraw = np.concatenate((bend[:,[0,1,2]], bend[:,[0,1,3]], np.stack(Ptri, axis=0)), axis=0)
    Triglraw = np.concatenate((bend[:,[0,1,2]], bend[:,[0,1,3]]), axis=0)
    Triglraw = np.sort(Triglraw, axis=1)
    Trigl = np.unique(Triglraw, axis=0)
    
#     Triglraw = np.vstack([bend[:, [0, 1, 2]], bend[:, [0, 1, 3]]])
#     Triglraw = np.sort(Triglraw, axis=1)
#     Trigl = np.unique(Triglraw, axis=0)


    # formulate connectivity matrix
    Comm = csr_matrix((Nn, Trigl.shape[0]))
    for i in range(Trigl.shape[0]):
        Comm[Trigl[i,:], i] = True
    
    # Comm = csr_matrix((np.ones(len(Trigl)), (Trigl[:,0], np.arange(len(Trigl)))), shape=(int(Nn+1), len(Trigl)))
    # Comm = Comm + csr_matrix((np.ones(len(Trigl)), (Trigl[:,1], np.arange(len(Trigl)))), shape=(int(Nn+1), len(Trigl)))
    # Comm = Comm + csr_matrix((np.ones(len(Trigl)), (Trigl[:,2], np.arange(len(Trigl)))), shape=(int(Nn+1), len(Trigl)))

    # search for fold lines
    Ge = Comm.T @ Comm
    mf, me = np.where(np.triu(Ge.toarray() == 2))

    # added this for ordering as in matlab code
    sorted_indices = np.argsort(me)
    me = me[sorted_indices]
    mf = mf[sorted_indices]


    fold = np.zeros((len(mf), 4), dtype=int)
    for i in range(len(mf)):
        link, ia, ib = np.intersect1d(Trigl[mf[i]], Trigl[me[i]], return_indices=True)
        oftpa = np.setdiff1d(np.arange(3), ia)
        oftpb = np.setdiff1d(np.arange(3), ib)
        fold[i] = [link[0], link[1], Trigl[mf[i], oftpa][0], Trigl[me[i], oftpb][0]]

    # from scipy.io import savemat
    # variables = {
    #     'fold2': fold,
    #     'link2': link,
    #     'ia2': ia,
    #     'ib2': ib,
    #     'oftpa2': oftpa,
    #     'oftpb2': oftpb,
    # }
    # # Save the variables to a .mat file
    # savemat('variables.mat', variables)

    # from scipy.io import savemat
    # variables = {
    #     'mf2': mf,
    #     'me2': me,
    # }
    # # Save the variables to a .mat file
    # savemat('variables.mat', variables)

    # fold = np.zeros((len(mf), 4))
    # for i in range(len(mf)):
    #     link, ia, ib = np.intersect1d(Trigl[mf[i]], Trigl[me[i]], return_indices=True)
    #     oftpa = [0, 1, 2][~np.isin([0, 1, 2], ia)]
    #     oftpb = [0, 1, 2][~np.isin([0, 1, 2], ib)]
    #     fold[i] = [link[0], link[1], Trigl[mf[i],oftpa], Trigl[me[i],oftpb]]
    
    fdandbd = np.sort(fold[:,:2], axis=1)
    onlybd = np.sort(bend[:,:2], axis=1)
    ibd = np.intersect1d(np.arange(len(fdandbd)), np.where(np.all(fdandbd == onlybd[:, None], axis=2))[1])
    fold = np.delete(fold, ibd, axis=0)
    
    # search for boundaries
    Edge = np.sort(np.concatenate((Trigl[:,:2], Trigl[:,1:], Trigl[:,[0,2]]), axis=0), axis=1)
    u, n, counts = np.unique(Edge, axis=0, return_inverse=True, return_counts=True)
    counts = np.bincount(n)
    bdry = u[counts == 1]    
    return fold, bdry, Trigl




# import numpy as np
# from scipy.sparse import lil_matrix, triu  

# def findfdbd(Panel, bend):
#     Nn = max(max(p) for p in Panel)+1   # sin +1 igual que matlab
#     # triangularization
#     Panelsize = [len(p) for p in Panel]
#     Ptri = [Panel[i] for i in range(len(Panel)) if Panelsize[i] == 3]
#     print(Ptri)
#     # Triglraw = np.vstack([bend[:, [0, 1, 2]], bend[:, [0, 1, 3]], Ptri])
#     Triglraw = np.vstack([bend[:, [0, 1, 2]], bend[:, [0, 1, 3]]])
#     Triglraw = np.sort(Triglraw, axis=1)
#     Trigl = np.unique(Triglraw, axis=0)    # -1 igual a matlab
#     # formulate connectivity matrix
#     Comm = lil_matrix((Nn, Trigl.shape[0]))
#     for i in range(Trigl.shape[0]):
#         indices = Trigl[i, :]
#         Comm[indices, i] = 1
    
#     # Comm = csr_matrix((np.ones(Trigl.shape[0]), (Trigl.ravel(), np.repeat(np.arange(Trigl.shape[0]), Trigl.shape[1]))), shape=(Nn, Trigl.shape[0]))
#     # search for fold lines
#     Ge = Comm.T @ Comm
#     mf, me = np.where(np.triu(Ge.toarray() == 2))
#     fold = np.zeros((len(mf), 4), dtype=int)
#     for i in range(len(mf)):
#         link, ia, ib = np.intersect1d(Trigl[mf[i]], Trigl[me[i]], return_indices=True)
#         oftpa = np.setdiff1d(np.arange(3), ia)
#         oftpb = np.setdiff1d(np.arange(3), ib)
#         fold[i] = [link[0], link[1], Trigl[mf[i], oftpa][0], Trigl[me[i], oftpb][0]]
#     fdandbd = np.sort(fold[:, :2], axis=1)
#     onlybd = np.sort(bend[:, :2], axis=1)
#     ibd, = np.where(np.isin(fdandbd, onlybd).all(axis=1))
#     fold = np.delete(fold, ibd, axis=0)

#     # search for boundaries
#     Edge = np.sort(np.vstack([Trigl[:, [0, 1]], Trigl[:, [1, 2]], Trigl[:, [2, 0]]]), axis=1)
#     u, n = np.unique(Edge, axis=0, return_counts=True)
#     bdry = u[n == 1]

#     return fold, bdry, Trigl
