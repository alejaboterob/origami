import numpy as np
from scipy import sparse

def GlobalK_fast_ver(Ui, Node, truss, angles):
    Nn = Node.shape[0]
    Nodenw = np.column_stack((
        Node[:, 0] + Ui[0::3],
        Node[:, 1] + Ui[1::3],
        Node[:, 2] + Ui[2::3]
    ))

    eDofb = (np.kron(truss['Bars'], 3*np.ones((1,3))) + np.tile([0,1,2], (truss['Bars'].shape[0], 2))).astype(int)
    du = Ui[eDofb[:, :3]] - Ui[eDofb[:, 3:6]]
    Ex = truss['B'] @ Ui / truss['L'] + 0.5 * np.sum(du**2, axis=1) / (truss['L']**2)
    
    Sx, Et, _ = truss['CM'](Ex)
    Duelem = np.column_stack((du, -du))
    Du = sparse.csr_matrix((Duelem.flatten(), (np.repeat(np.arange(len(Et)), 6), eDofb.flatten())), shape=(len(Et), len(Ui)))
    Fx = Sx * truss['A']
    IFb = (np.sum(truss['B'].T * Fx[:, np.newaxis], axis=1) + np.sum(Du.T.multiply(Fx / truss['L']), axis=1).T)
    
    Kel = (truss['B'].T @ sparse.diags(Et * truss['A'] / truss['L']) @ truss['B']).T
    K1 = Du.T @ sparse.diags(Et * truss['A'] / (truss['L']**2)) @ truss['B'] + truss['B'].T @ sparse.diags(Et * truss['A'] / (truss['L']**2)) @ Du
    K2 = Du.T @ sparse.diags(Et * truss['A'] / (truss['L']**3)) @ Du

    G = sparse.csr_matrix((np.concatenate([np.ones(len(Et)), -np.ones(len(Et))]), 
                        (np.concatenate([np.arange(len(Et)), np.arange(len(Et))]), truss['Bars'].T.flatten())), 
                        shape=(len(Et), Nn))
    
    ia, ja, sa = sparse.find(G.T @ sparse.diags(Fx / truss['L']) @ G)
    ik = np.outer(3*(ia), np.ones(3)) + np.arange(0, 3)
    jk = np.outer(3*(ja), np.ones(3)) + np.arange(0, 3)
    Kg = sparse.csr_matrix((np.tile(sa, (3,1)).T.flatten(), (jk.flatten(), ik.flatten())), shape=(3*Nn, 3*Nn))    
    Kg = 0.5 * (Kg + Kg.T)
    Kb = (Kel + K1 + K2) + Kg

    rotspr = np.vstack((angles['bend'], angles['fold']))
    
    if angles['pb0'].size == 0:
        h0 = angles['pf0']
    else:
        h0 = np.hstack((angles['pb0'], angles['pf0']))
    
    if angles['kpb'].size == 0:
        kpi = angles['kpf']
    else:
        kpi = np.hstack((angles['kpb'], angles['kpf']))

    eDofd = np.kron(rotspr, 3*np.ones((1,3))) + np.tile([0,1,2], (rotspr.shape[0], 4))
    rkj = (Nodenw[rotspr[:, 1]] - Nodenw[rotspr[:, 0]])
    rij = (Nodenw[rotspr[:, 2]] - Nodenw[rotspr[:, 0]])
    rkl = (Nodenw[rotspr[:, 1]] - Nodenw[rotspr[:, 3]])
    

    rmj = np.cross(rij, rkj)
    rnk = np.cross(rkj, rkl)
    dt_rnkrij = np.sum(rnk * rij, axis=1)
    sgn = np.where(np.abs(dt_rnkrij) > 1e-8, np.sign(dt_rnkrij), 1)
    dt_rmjrnk = np.sum(rmj * rnk, axis=1)
    rmj2 = np.sum(rmj**2, axis=1)
    norm_rmj = np.sqrt(rmj2)
    rkj2 = np.sum(rkj**2, axis=1)
    norm_rkj = np.sqrt(rkj2)
    rnk2 = np.sum(rnk**2, axis=1)
    norm_rnk = np.sqrt(rnk2)
    he = np.array([sgn * np.real(np.arccos(dt_rmjrnk / (norm_rmj * norm_rnk)))])
    he[he < 0] += 2*np.pi
    Rspr, Kspr, _ = angles['CM'](he, h0, kpi, truss['L'][:rotspr.shape[0]])

    dt_rijrkj = np.sum(rij * rkj, axis=1)
    dt_rklrkj = np.sum(rkl * rkj, axis=1)

    di = np.einsum('mi, m->mi', rmj, (norm_rkj / rmj2))
    dl = np.einsum('mi, m->mi', -rnk, (norm_rkj / rnk2))
    dj = np.einsum('mi, m->mi', di, (dt_rijrkj / rkj2 - 1)) - np.einsum('mi, m->mi', dl, (dt_rklrkj / rkj2))
    dk = np.einsum('mi, m->mi', -di, (dt_rijrkj / rkj2)) + np.einsum('mi, m->mi', dl, (dt_rklrkj / rkj2 - 1))
    Jhe_dense = np.hstack((dj, dk, di, dl))
    Jhe = sparse.csr_matrix((Jhe_dense.flatten(), (eDofd.flatten(), np.repeat(np.arange(len(he.T)), 12))), shape=(len(Ui), len(he.T)))
    IFbf = np.sum(Jhe.multiply(Rspr), axis=1)
    IF = IFb.T + IFbf

    # dii = -  np.multiply( np.multiply(np.transpose(rmj[:, np.newaxis, :], (0, 2, 1)), np.transpose(np.cross(rkj, rmj)[:, np.newaxis, :], (2, 0, 1)))   + \
    #             np.transpose(np.multiply(np.transpose(rmj[:, np.newaxis, :], (0, 2, 1)), np.transpose(np.cross(rkj, rmj)[:, np.newaxis, :], (2, 0, 1))), (1, 0, 2)), 
    #             np.transpose(np.repeat( norm_rkj / (rmj2**2), 3)[:, np.newaxis, np.newaxis], (2, 0, 1)))

    # dii = -bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj,rmj),[3 1 2]))+...
    #                      permute(bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj,rmj),[3 1 2])),[2 1 3]),...
    #               permute(norm_rkj./(rmj2.^2),[3 1 2]));

#     dii = - (np.einsum('mij,kim->mjk', rmj[:, np.newaxis, :], np.cross(rkj, rmj)[:, np.newaxis, :]) + \
#             np.einsum('mij,kim->jmk', rmj[:, np.newaxis, :], np.cross(rkj, rmj)[:, np.newaxis, :])) * \
#             (norm_rkj / (rmj2**2))
    
    dii = - (np.einsum('mi,mj->mij', rmj, np.cross(rkj, rmj)) + \
         np.einsum('mi,mj->mji', rmj, np.cross(rkj, rmj))) * \
        (norm_rkj[:, np.newaxis, np.newaxis] / (rmj2[:, np.newaxis, np.newaxis] ** 2))
    
    dij = - np.einsum('mi,mj->mij', rmj, rkj) * (1 / (rmj2[:, np.newaxis, np.newaxis] * norm_rkj[:, np.newaxis, np.newaxis])) + \
                (np.einsum('mi,mj->mij', rmj, np.cross(rkj - rij, rmj)) + \
                np.einsum('mi,mj->mji', rmj, np.cross(rkj - rij, rmj))) * \
                (norm_rkj[:, np.newaxis, np.newaxis] / (rmj2[:, np.newaxis, np.newaxis]**2))


#     dij = - np.einsum('mij,kim->mjk', rmj[:, np.newaxis, :], rkj[:, np.newaxis, :]) * (1 / (rmj2 * norm_rkj)) + \
#             (np.einsum('mij,kim->mjk', rmj[:, np.newaxis, :], np.cross(rkj - rij, rmj)[:, np.newaxis, :]) + \
#             np.einsum('mij,kim->jmk', rmj[:, np.newaxis, :], np.cross(rkj - rij, rmj)[:, np.newaxis, :])) * \
#             (norm_rkj / (rmj2**2))

    # dij = -bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(rkj,[3 1 2])),...
    #               permute(1./(rmj2.*norm_rkj),[3 1 2]))...
    #       +bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj-rij,rmj),[3 1 2]))+...
    #                             permute(bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj-rij,rmj),[3 1 2])),[2 1 3]),...
    #               permute(norm_rkj./(rmj2.^2),[3 1 2]));

    # dij = - np.multiply( np.multiply(np.transpose(rmj[:, np.newaxis, :], (0, 2, 1)), np.transpose(rkj[:, np.newaxis, :], (2, 0, 1))), 
    #                     np.transpose (np.repeat(1 / (rmj2 * norm_rkj), 3)[:, np.newaxis, np.newaxis], (2, 0, 1))) + \
    #         np.multiply( np.multiply(np.transpose(rmj[:, np.newaxis, :], (0, 2, 1)), np.transpose(np.cross(rkj - rij, rmj)[:, np.newaxis, :], (2, 0, 1))) + \
    #                     np.transpose(np.multiply(np.transpose(rmj[:, np.newaxis, :], (0, 2, 1)), np.transpose(np.cross(rkj - rij, rmj)[:, np.newaxis, :], (2, 0, 1))), (1, 0, 2)),  
    #                     np.transpose(np.repeat(norm_rkj / (rmj2**2), 3)[:, np.newaxis, np.newaxis], (2, 0, 1)))

    dik = np.einsum('mi,mj->mij', rmj, rkj) * (1 / (rmj2[:, np.newaxis, np.newaxis] * norm_rkj[:, np.newaxis, np.newaxis])) + \
                (np.einsum('mi,mj->mij', rmj, np.cross(rij, rmj)) + \
                np.einsum('mi,mj->mji', rmj, np.cross(rij, rmj))) * \
                (norm_rkj[:, np.newaxis, np.newaxis] / (rmj2[:, np.newaxis, np.newaxis]**2))
    

#     dik = np.einsum('mij,kim->mjk', rmj[:, np.newaxis, :], rkj[:, np.newaxis, :]) * (1 / (rmj2 * norm_rkj)) + \
#             (np.einsum('mij,kim->mjk', rmj[:, np.newaxis, :], np.cross(rij, rmj)[:, np.newaxis, :]) + \
#             np.einsum('mij,kim->jmk', rmj[:, np.newaxis, :], np.cross(rij, rmj)[:, np.newaxis, :])) * \
#             (norm_rkj / (rmj2**2))

    #   dik =  bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(rkj,[3 1 2])),...
    #                   permute(1./(rmj2.*norm_rkj),[3 1 2]))...
    #           +bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rij,rmj),[3 1 2]))+...
    #                              permute(bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rij,rmj),[3 1 2])),[2 1 3]),...
    #                       permute(norm_rkj./(rmj2.^2),[3 1 2]));
    dll = (np.einsum('mi,mj->mij', rnk, np.cross(rkj, rnk)) + \
           np.einsum('mi,mj->mji', rnk, np.cross(rkj, rnk))) * \
          (norm_rkj[:, np.newaxis, np.newaxis] / (rnk2[:, np.newaxis, np.newaxis]**2))
    

#     dll = (np.einsum('ijk,jkl->kli', rnk[:, np.newaxis, :], np.cross(rkj, rnk)[:, np.newaxis, :]) +
#             np.einsum('ijk,jkl->lki', rnk[:, np.newaxis, :], np.cross(rkj, rnk)[:, np.newaxis, :])) * \
#             (norm_rkj / (rnk2**2))
                

    #     dll = bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj,rnk),[3 1 2]))+...
    #                         permute(bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj,rnk),[3 1 2])),[2 1 3]),...
    #                  permute(norm_rkj./(rnk2.^2),[3 1 2]));

    dlk = - np.einsum('mi,mj->mij', rnk, rkj) * (1 / (rnk2[:, np.newaxis, np.newaxis] * norm_rkj[:, np.newaxis, np.newaxis])) - \
                (np.einsum('mi,mj->mij', rnk, np.cross(rkj - rkl, rnk)) + \
                np.einsum('mi,mj->mji', rnk, np.cross(rkj - rkl, rnk))) * \
                (norm_rkj[:, np.newaxis, np.newaxis] / (rnk2[:, np.newaxis, np.newaxis]**2))

#     dlk = - np.einsum('ijk,jkl->lki', rnk[:, np.newaxis, :], rkj[:, np.newaxis, :]) * (1 / (rnk2 * norm_rkj))  \
#             - (np.einsum('ijk,jkl->lki', rnk[:, np.newaxis, :], np.cross(rkj - rkl, rnk)[:, np.newaxis, :]) +
#             np.einsum('ijk,jkl->kli', rnk[:, np.newaxis, :], np.cross(rkj - rkl, rnk)[:, np.newaxis, :])) * \
#             (norm_rkj / (rnk2**2))
    
    #     dlk = -bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(rkj,[3 1 2])),...
    #                   permute(1./(rnk2.*norm_rkj),[3 1 2]))...
    #           -bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj-rkl,rnk),[3 1 2]))+...
    #                          permute(bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj-rkl,rnk),[3 1 2])),[2 1 3]),...
    #                   permute(norm_rkj./(rnk2.^2),[3 1 2]));

    dlj = np.einsum('mi,mj->mij', rnk, rkj) * (1 / (rnk2[:, np.newaxis, np.newaxis] * norm_rkj[:, np.newaxis, np.newaxis])) - \
                (np.einsum('mi,mj->mij', rnk, np.cross(rkl, rnk)) + \
                np.einsum('mi,mj->mji', rnk, np.cross(rkl, rnk))) * \
                (norm_rkj[:, np.newaxis, np.newaxis] / (rnk2[:, np.newaxis, np.newaxis]**2))
    

#     dlj = np.einsum('ijk,jkl->lki', rnk[:, np.newaxis, :], rkj[:, np.newaxis, :]) * (1 / (rnk2 * norm_rkj)) - \
#             (np.einsum('ijk,jkl->lki', rnk[:, np.newaxis, :], np.cross(rkl, rnk)[:, np.newaxis, :]) + \
#             np.einsum('ijk,jkl->kli', rnk[:, np.newaxis, :], np.cross(rkl, rnk)[:, np.newaxis, :])) * \
#             (norm_rkj / (rnk2**2))

    #     dlj =  bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(rkj,[3 1 2])),permute(1./(rnk2.*norm_rkj),[3 1 2]))...
    #           -bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkl,rnk),[3 1 2]))+...
    #                          permute(bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkl,rnk),[3 1 2])),[2 1 3]),...
    #                   permute(norm_rkj./(rnk2.^2),[3 1 2]));

    dT1jj = np.einsum('mi,m->mi', np.einsum('mi,m->mi', rkj, (-1 + 2 * dt_rijrkj / rkj2)) - rij, 1/ rkj2) 
    dT2jj = np.einsum('mi,m->mi', np.einsum('mi,m->mi', rkj, (2 * dt_rklrkj / rkj2)) - rkl, 1/ rkj2)
    djj = np.einsum('mi, mj->mij', di, dT1jj) + \
        np.einsum('mij, m->mij', dij, dt_rijrkj / rkj2 - 1) - \
        np.einsum('mi, mj->mij', dl, dT2jj) - \
        np.einsum('mij, m->mij', dlj, dt_rklrkj / rkj2)


#     dT1jj = ((rkj * (-1 + 2 * dt_rijrkj / rkj2) - rij) / rkj2  ) [:, :, np.newaxis]
#     dT2jj = ((rkj * (2 * dt_rklrkj / rkj2) - rkl) / rkj2)[: , :, np.newaxis]


#     djj = np.einsum('mi,ijk->jik', di, dT1jj) + \
#             np.einsum('mij,jik->mik', dij, np.array([dt_rijrkj / rkj2 - 1])[:, np.newaxis, np.newaxis])  - \
#             np.einsum('mi,ijk->jik', dl, dT2jj) - \
#             np.einsum('mij,jik->mik', dlj, np.array([dt_rklrkj / rkj2])[:, np.newaxis, np.newaxis])
    
    #     dT1jj = bsxfun(@times,bsxfun(@times,rkj,(-1+2*dt_rijrkj./rkj2))-rij,1./(rkj2));
    #     dT2jj = bsxfun(@times,bsxfun(@times,rkj,(2*dt_rklrkj./rkj2))-rkl,1./(rkj2));
    #     djj = bsxfun(@times,permute(di,[1 3 2]),permute(dT1jj,[3 1 2]))+...
    #           bsxfun(@times,dij,permute((dt_rijrkj./rkj2)-1,[3 1 2]))-...
    #           bsxfun(@times,permute(dl,[1 3 2]),permute(dT2jj,[3 1 2]))-...
    #           bsxfun(@times,dlj,permute(dt_rklrkj./rkj2,[3 1 2]));

    dT1jk = np.einsum('mi,m->mi', np.einsum('mi,m->mi', rkj, (-2 * dt_rijrkj / rkj2)) + rij, 1/ rkj2)
    dT2jk = np.einsum('mi,m->mi', np.einsum('mi,m->mi', rkj, (1 - 2 * dt_rklrkj / rkj2)) + rkl, 1/ rkj2)
    djk = np.einsum('mi, mj->mij', di, dT1jk) + \
        np.einsum('mij, m->mij', dik, dt_rijrkj / rkj2 - 1) - \
        np.einsum('mi, mj->mij', dl, dT2jk) - \
        np.einsum('mij, m->mij', dlk, dt_rklrkj / rkj2)
    

#     dT1jk = ((rkj * (-2 * dt_rijrkj / rkj2) + rij) / rkj2)[:, :, np.newaxis]
#     dT2jk = ((rkj * (1 - 2 * dt_rklrkj / rkj2) + rkl) / rkj2)[:, :, np.newaxis]
#     djk = np.einsum('mi,ijk->jik', di, dT1jk) + \
#             np.einsum('mij,jik->mik', dik, np.array([dt_rijrkj / rkj2 - 1])[:, np.newaxis, np.newaxis])  - \
#             np.einsum('mi,ijk->jik', dl, dT2jk) - \
#             np.einsum('mij,jik->mik', dlk, np.array([dt_rklrkj / rkj2])[:, np.newaxis, np.newaxis])

    #     dT1jk = bsxfun(@times,bsxfun(@times,rkj,(-2*dt_rijrkj./rkj2))+rij,1./(rkj2));
    #     dT2jk = bsxfun(@times,bsxfun(@times,rkj,(1-2*dt_rklrkj./rkj2))+rkl,1./(rkj2));
    #     djk = bsxfun(@times,permute(di,[1 3 2]),permute(dT1jk,[3 1 2]))+...
    #           bsxfun(@times,dik,permute((dt_rijrkj./rkj2)-1,[3 1 2]))-...
    #           bsxfun(@times,permute(dl,[1 3 2]),permute(dT2jk,[3 1 2]))-...
    #           bsxfun(@times,dlk,permute(dt_rklrkj./rkj2,[3 1 2]));

    dT1kk = dT2jk
    dT2kk = dT1jk
    dkk = np.einsum('mi, mj->mij', dl, dT1kk) + \
        np.einsum('mij, m->mij', dlk, dt_rklrkj / rkj2 - 1) - \
        np.einsum('mi, mj->mij', di, dT2kk) - \
        np.einsum('mij, m->mij', dik, dt_rijrkj / rkj2)
    


#     dkk = np.einsum('mi,ijk->jik', dl, dT1kk) + \
#             np.einsum('mij,jik->mik', dlk, np.array([dt_rklrkj / rkj2 - 1])[:, np.newaxis, np.newaxis])  - \
#             np.einsum('mi,ijk->jik', di, dT2kk) - \
#             np.einsum('mij,jik->mik', dik, np.array([dt_rijrkj / rkj2])[:, np.newaxis, np.newaxis])
    

    #     dT1kk = dT2jk; dT2kk = dT1jk;
    #     dkk = bsxfun(@times,permute(dl,[1 3 2]),permute(dT1kk,[3 1 2]))+...
    #           bsxfun(@times,dlk,permute(dt_rklrkj./rkj2-1,[3 1 2]))-...
    #           bsxfun(@times,permute(di,[1 3 2]),permute(dT2kk,[3 1 2]))-...
    #           bsxfun(@times,dik,permute(dt_rijrkj./rkj2,[3 1 2]));
    
    Hp = np.zeros((12, 12, len(he.T)))
    Hp[0:3, 0:3, :] = djj.T
    Hp[6:9, 6:9, :] = dii.T
    Hp[3:6, 3:6, :] = dkk.T
    Hp[9:12, 9:12, :] = dll.T
    Hp[0:3, 3:6, :] = np.transpose(djk.T, (1, 0, 2))
    Hp[3:6, 0:3, :] = djk.T
    Hp[6:9, 0:3, :] = np.transpose(dij.T, (1, 0, 2))
    Hp[0:3, 6:9, :] = dij.T
    Hp[9:12, 0:3, :] = np.transpose(dlj.T, (1, 0, 2))
    Hp[0:3, 9:12, :] = dlj.T
    Hp[6:9, 3:6, :] = np.transpose(dik.T, (1, 0, 2))
    Hp[3:6, 6:9, :] = dik.T
    Hp[9:12, 3:6, :] = np.transpose(dlk.T, (1, 0, 2))
    Hp[3:6, 9:12, :] = dlk.T

    Khe_dense = np.einsum('ijk,j->ijk', 
                        np.einsum('ij,jk->jki', Jhe_dense, Jhe_dense),
                        Kspr) + \
                np.einsum('ijk,j->ijk', Hp, Rspr)
    
    dof_ind = np.tile(eDofd, (12,1))
    Kbf = sparse.coo_matrix((Khe_dense.flatten(), 
                            (dof_ind.T.flatten(), dof_ind.flatten())), 
                            shape=(3*Nn, 3*Nn))

    K = Kb + Kbf
    K = (K + K.T) / 2

    return IF, K
