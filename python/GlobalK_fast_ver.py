import numpy as np
from scipy import sparse

# Assuming icross is a custom function for cross product, you might need to implement it
def icross(a, b):
    return np.cross(a, b, axisa=0, axisb=0, axisc=0)

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
    
    Kel = truss['B'].T @ sparse.diags(Et * truss['A'] / truss['L']) @ truss['B']
    K1 = Du.T @ sparse.diags(Et * truss['A'] / (truss['L']**2)) @ truss['B'] + truss['B'].T @ sparse.diags(Et * truss['A'] / (truss['L']**2)) @ Du
    K2 = Du.T @ sparse.diags(Et * truss['A'] / (truss['L']**3)) @ Du
    G = sparse.csr_matrix((np.concatenate([np.ones(len(Et)), -np.ones(len(Et))]), 
                           (np.concatenate([np.arange(len(Et)), np.arange(len(Et))]), truss['Bars'].flatten())), 
                          shape=(len(Et), Nn))
    
    ia, ja, sa = sparse.find(G.T @ sparse.diags(Fx / truss['L']) @ G)
    ik = np.outer(3*(ia-1), np.ones(3)) + np.arange(1, 4)
    jk = np.outer(3*(ja-1), np.ones(3)) + np.arange(1, 4)
    Kg = sparse.csr_matrix((np.tile(sa, (3, 1)).flatten(), (ik.flatten(), jk.flatten())), shape=(3*Nn, 3*Nn))
    Kg = 0.5 * (Kg + Kg.T)
    Kb = (Kel + K1 + K2) + Kg

    rotspr = np.vstack((angles['bend'], angles['fold']))
    
    if angles['pb0'].size == 0:
        h0 = angles['pf0']
    else:
        h0 = np.vstack((angles['pb0'], angles['pf0']))
    
    if angles['kpb'].size == 0:
        kpi = angles['kpf']
    else:
        kpi = np.column_stack((angles['kpb'], angles['kpf']))

    eDofd = np.kron(rotspr, 3*np.ones((1,3))) + np.tile([0,1,2], (rotspr.shape[0], 4))
    rkj = (Nodenw[rotspr[:, 1]] - Nodenw[rotspr[:, 0]]).T
    rij = (Nodenw[rotspr[:, 2]] - Nodenw[rotspr[:, 0]]).T
    rkl = (Nodenw[rotspr[:, 1]] - Nodenw[rotspr[:, 3]]).T
    
    def icross(a, b):
        return np.cross(a.T, b.T).T
    
    rmj = icross(rij, rkj)
    rnk = icross(rkj, rkl)
    dt_rnkrij = np.sum(rnk * rij, axis=0)
    sgn = np.where(np.abs(dt_rnkrij) > 1e-8, np.sign(dt_rnkrij), 1)
    dt_rmjrnk = np.sum(rmj * rnk, axis=0)
    rmj2 = np.sum(rmj**2, axis=0)
    norm_rmj = np.sqrt(rmj2)
    rkj2 = np.sum(rkj**2, axis=0)
    norm_rkj = np.sqrt(rkj2)
    rnk2 = np.sum(rnk**2, axis=0)
    norm_rnk = np.sqrt(rnk2)
    he = sgn * np.real(np.arccos(dt_rmjrnk / (norm_rmj * norm_rnk)))
    he[he < 0] += 2*np.pi
    Rspr, Kspr, _ = angles['CM'](he, h0, kpi, truss['L'][:rotspr.shape[0]])

    dt_rijrkj = np.sum(rij * rkj, axis=0)
    dt_rklrkj = np.sum(rkl * rkj, axis=0)

    di = rmj * (norm_rkj / rmj2)
    dl = -rnk * (norm_rkj / rnk2)
    dj = di * (dt_rijrkj / rkj2 - 1) - dl * (dt_rklrkj / rkj2)
    dk = -di * (dt_rijrkj / rkj2) + dl * (dt_rklrkj / rkj2 - 1)
    Jhe_dense = np.vstack((dj, dk, di, dl))
    Jhe = sparse.csr_matrix((Jhe_dense.flatten(), (eDofd.flatten(), np.repeat(np.arange(len(he)), 12))), shape=(len(Ui), len(he)))
    IFbf = np.sum(Jhe.multiply(Rspr), axis=1)
    IF = IFb.T + IFbf

    dii = -np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(icross(rkj, rmj), (2, 0, 1))) +
                 np.transpose(np.einsum('ijk,kjl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(icross(rkj, rmj), (2, 0, 1))), (1, 0, 2)),
                 np.transpose(norm_rkj / (rmj2 ** 2), (2, 0, 1)))

    dij = -np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(rkj, (2, 0, 1))),
                    np.transpose(1. / (rmj2 * norm_rkj), (2, 0, 1))) + \
        np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(icross(rkj - rij, rmj), (2, 0, 1))) +
                    np.transpose(np.einsum('ijk,kjl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(icross(rkj - rij, rmj), (2, 0, 1))), (1, 0, 2)),
                    np.transpose(norm_rkj / (rmj2 ** 2), (2, 0, 1)))

    dik = np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(rkj, (2, 0, 1))),
                    np.transpose(1. / (rmj2 * norm_rkj), (2, 0, 1))) + \
        np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(icross(rij, rmj), (2, 0, 1))) +
                    np.transpose(np.einsum('ijk,kjl->ijl', np.transpose(rmj, (0, 2, 1)), np.transpose(icross(rij, rmj), (2, 0, 1))), (1, 0, 2)),
                    np.transpose(norm_rkj / (rmj2 ** 2), (2, 0, 1)))

    dll = np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(icross(rkj, rnk), (2, 0, 1))) +
                    np.transpose(np.einsum('ijk,kjl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(icross(rkj, rnk), (2, 0, 1))), (1, 0, 2)),
                    np.transpose(norm_rkj / (rnk2 ** 2), (2, 0, 1)))

    dlk = -np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(rkj, (2, 0, 1))),
                    np.transpose(1. / (rnk2 * norm_rkj), (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(icross(rkj - rkl, rnk), (2, 0, 1))) +
                    np.transpose(np.einsum('ijk,kjl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(icross(rkj - rkl, rnk), (2, 0, 1))), (1, 0, 2)),
                    np.transpose(norm_rkj / (rnk2 ** 2), (2, 0, 1)))

    dlj = np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(rkj, (2, 0, 1))),
                    np.transpose(1. / (rnk2 * norm_rkj), (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', np.einsum('ijk,jkl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(icross(rkl, rnk), (2, 0, 1))) +
                    np.transpose(np.einsum('ijk,kjl->ijl', np.transpose(rnk, (0, 2, 1)), np.transpose(icross(rkl, rnk), (2, 0, 1))), (1, 0, 2)),
                    np.transpose(norm_rkj / (rnk2 ** 2), (2, 0, 1)))

    dT1jj = np.einsum('ijk,jk->ij', rkj, (-1 + 2 * dt_rijrkj / rkj2)) - rij
    dT2jj = np.einsum('ijk,jk->ij', rkj, (2 * dt_rklrkj / rkj2)) - rkl
    djj = np.einsum('ijk,kjl->ijl', np.transpose(di, (0, 2, 1)), np.transpose(dT1jj, (2, 0, 1))) + \
        np.einsum('ijk,kjl->ijl', dij, np.transpose((dt_rijrkj / rkj2) - 1, (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', np.transpose(dl, (0, 2, 1)), np.transpose(dT2jj, (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', dlj, np.transpose(dt_rklrkj / rkj2, (2, 0, 1)))

    dT1jk = np.einsum('ijk,jk->ij', rkj, (-2 * dt_rijrkj / rkj2)) + rij
    dT2jk = np.einsum('ijk,jk->ij', rkj, (1 - 2 * dt_rklrkj / rkj2)) + rkl
    djk = np.einsum('ijk,kjl->ijl', np.transpose(di, (0, 2, 1)), np.transpose(dT1jk, (2, 0, 1))) + \
        np.einsum('ijk,kjl->ijl', dik, np.transpose((dt_rijrkj / rkj2) - 1, (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', np.transpose(dl, (0, 2, 1)), np.transpose(dT2jk, (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', dlk, np.transpose(dt_rklrkj / rkj2, (2, 0, 1)))

    dT1kk = dT2jk
    dT2kk = dT1jk
    dkk = np.einsum('ijk,kjl->ijl', np.transpose(dl, (0, 2, 1)), np.transpose(dT1kk, (2, 0, 1))) + \
        np.einsum('ijk,kjl->ijl', dlk, np.transpose(dt_rklrkj / rkj2 - 1, (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', np.transpose(di, (0, 2, 1)), np.transpose(dT2kk, (2, 0, 1))) - \
        np.einsum('ijk,kjl->ijl', dik, np.transpose(dt_rijrkj / rkj2, (2, 0, 1)))




    dii = -np.einsum('ijk,jlk->ilk', 
                 np.einsum('ijk,jlk->ilk', rmj[:, np.newaxis, :], icross(rkj, rmj).T[:, :, np.newaxis]),
                 (norm_rkj / (rmj2**2))[np.newaxis, :, :])

    dij = -np.einsum('ijk,jlk->ilk',
                    np.einsum('ijk,jlk->ilk', rmj[:, np.newaxis, :], rkj.T[:, :, np.newaxis]),
                    (1 / (rmj2 * norm_rkj))[np.newaxis, :, :]) + \
        np.einsum('ijk,jlk->ilk',
                np.einsum('ijk,jlk->ilk', rmj[:, np.newaxis, :], icross(rkj - rij, rmj).T[:, :, np.newaxis]),
                (norm_rkj / (rmj2**2))[np.newaxis, :, :])

    dik = np.einsum('ijk,jlk->ilk',
                    np.einsum('ijk,jlk->ilk', rmj[:, np.newaxis, :], rkj.T[:, :, np.newaxis]),
                    (1 / (rmj2 * norm_rkj))[np.newaxis, :, :]) + \
        np.einsum('ijk,jlk->ilk',
                np.einsum('ijk,jlk->ilk', rmj[:, np.newaxis, :], icross(rij, rmj).T[:, :, np.newaxis]),
                (norm_rkj / (rmj2**2))[np.newaxis, :, :])

    # dil = np.zeros((3, 3, len(he)))

    dll = np.einsum('ijk,jlk->ilk',
                    np.einsum('ijk,jlk->ilk', rnk[:, np.newaxis, :], icross(rkj, rnk).T[:, :, np.newaxis]),
                    (norm_rkj / (rnk2**2))[np.newaxis, :, :])

    dlk = -np.einsum('ijk,jlk->ilk',
                    np.einsum('ijk,jlk->ilk', rnk[:, np.newaxis, :], rkj.T[:, :, np.newaxis]),
                    (1 / (rnk2 * norm_rkj))[np.newaxis, :, :]) - \
        np.einsum('ijk,jlk->ilk',
                np.einsum('ijk,jlk->ilk', rnk[:, np.newaxis, :], icross(rkj - rkl, rnk).T[:, :, np.newaxis]),
                (norm_rkj / (rnk2**2))[np.newaxis, :, :])

    dlj = np.einsum('ijk,jlk->ilk',
                    np.einsum('ijk,jlk->ilk', rnk[:, np.newaxis, :], rkj.T[:, :, np.newaxis]),
                    (1 / (rnk2 * norm_rkj))[np.newaxis, :, :]) - \
        np.einsum('ijk,jlk->ilk',
                np.einsum('ijk,jlk->ilk', rnk[:, np.newaxis, :], icross(rkl, rnk).T[:, :, np.newaxis]),
                (norm_rkj / (rnk2**2))[np.newaxis, :, :])

    dT1jj = np.einsum('ijk,jk->ijk', rkj * (-1 + 2 * dt_rijrkj / rkj2) - rij, 1 / rkj2)
    dT2jj = np.einsum('ijk,jk->ijk', rkj * (2 * dt_rklrkj / rkj2) - rkl, 1 / rkj2)
    djj = np.einsum('ijk,jlk->ilk', di[:, np.newaxis, :], dT1jj.transpose(2, 0, 1)[:, :, np.newaxis]) + \
        np.einsum('ijk,jk->ijk', dij, (dt_rijrkj / rkj2 - 1)[np.newaxis, :]) - \
        np.einsum('ijk,jlk->ilk', dl[:, np.newaxis, :], dT2jj.transpose(2, 0, 1)[:, :, np.newaxis]) - \
        np.einsum('ijk,jk->ijk', dlj, (dt_rklrkj / rkj2)[np.newaxis, :])

    dT1jk = np.einsum('ijk,jk->ijk', rkj * (-2 * dt_rijrkj / rkj2) + rij, 1 / rkj2)
    dT2jk = np.einsum('ijk,jk->ijk', rkj * (1 - 2 * dt_rklrkj / rkj2) + rkl, 1 / rkj2)
    djk = np.einsum('ijk,jlk->ilk', di[:, np.newaxis, :], dT1jk.transpose(2, 0, 1)[:, :, np.newaxis]) + \
        np.einsum('ijk,jk->ijk', dik, (dt_rijrkj / rkj2 - 1)[np.newaxis, :]) - \
        np.einsum('ijk,jlk->ilk', dl[:, np.newaxis, :], dT2jk.transpose(2, 0, 1)[:, :, np.newaxis]) - \
        np.einsum('ijk,jk->ijk', dlk, (dt_rklrkj / rkj2)[np.newaxis, :])

    dT1kk = dT2jk
    dT2kk = dT1jk
    dkk = np.einsum('ijk,jlk->ilk', dl[:, np.newaxis, :], dT1kk.transpose(2, 0, 1)[:, :, np.newaxis]) + \
        np.einsum('ijk,jk->ijk', dlk, (dt_rklrkj / rkj2 - 1)[np.newaxis, :]) - \
        np.einsum('ijk,jlk->ilk', di[:, np.newaxis, :], dT2kk.transpose(2, 0, 1)[:, :, np.newaxis]) - \
        np.einsum('ijk,jk->ijk', dik, (dt_rijrkj / rkj2)[np.newaxis, :])

    Hp = np.zeros((12, 12, len(he)))
    Hp[0:3, 0:3, :] = djj
    Hp[6:9, 6:9, :] = dii
    Hp[3:6, 3:6, :] = dkk
    Hp[9:12, 9:12, :] = dll
    Hp[0:3, 3:6, :] = djk
    Hp[3:6, 0:3, :] = np.transpose(djk, (1, 0, 2))
    Hp[6:9, 0:3, :] = dij
    Hp[0:3, 6:9, :] = np.transpose(dij, (1, 0, 2))
    Hp[9:12, 0:3, :] = dlj
    Hp[0:3, 9:12, :] = np.transpose(dlj, (1, 0, 2))
    Hp[6:9, 3:6, :] = dik
    Hp[3:6, 6:9, :] = np.transpose(dik, (1, 0, 2))
    Hp[9:12, 3:6, :] = dlk
    Hp[3:6, 9:12, :] = np.transpose(dlk, (1, 0, 2))

    Khe_dense = np.einsum('ijk,jlk->ilk', 
                        np.einsum('ijk,jlk->ilk', Jhe_dense[:, np.newaxis, :], Jhe_dense.T[:, :, np.newaxis]),
                        Kspr[np.newaxis, :, :]) + \
                np.einsum('ijk,jk->ijk', Hp, Rspr[np.newaxis, :])

    dof_ind1 = np.tile(eDofd, (1, 1, 12)).transpose(0, 2, 1)
    dof_ind2 = dof_ind1.transpose(1, 0, 2)
    Kbf = sparse.coo_matrix((Khe_dense.flatten(), 
                            (dof_ind1.flatten(), dof_ind2.flatten())), 
                            shape=(3*Nn, 3*Nn))

    K = Kb + Kbf
    K = (K + K.T) / 2

    return IF, Kb
