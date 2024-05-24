import numpy as np

def icross(a, b):
    return np.cross(a, b, axis=0)

def FoldKe(Cood, List, kpi, h0=None, L0=None, CM=None):
    rkj = (Cood[List[1]] - Cood[List[0]]).T
    rij = (Cood[List[2]] - Cood[List[0]]).T
    rkl = (Cood[List[1]] - Cood[List[3]]).T
    rmj = icross(rij, rkj)
    rnk = icross(rkj, rkl)
    sgn = ((abs(np.dot(rnk.T, rij)) > 1e-8) * np.sign(np.dot(rnk.T, rij)) + 
           (abs(np.dot(rnk.T, rij)) <= 1e-8) * 1)
    he = np.real(np.arccos(np.dot(rmj.T, rnk) / (np.linalg.norm(rmj) * np.linalg.norm(rnk))))
    he = np.real(sgn * he)
    if he < 0:
        he = 2 * np.pi + he

    Rhe = None
    Khe = None
    if h0 is not None and L0 is not None and CM is not None:
        Rspr, Kspr = CM(he, h0, kpi, L0)
    
        di = np.linalg.norm(rkj) / np.dot(rmj.T, rmj) * rmj
        dl = -np.linalg.norm(rkj) / np.dot(rnk.T, rnk) * rnk
        dj = (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * di - np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dl
        dk = -np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * di + (np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) - 1) * dl

        Jhe = np.vstack([dj, dk, di, dl])
        Rhe = (Rspr * Jhe).flatten()

        dii = -np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (np.outer(rmj, np.cross(rkj, rmj).T) + np.outer(rmj, np.cross(rkj, rmj).T).T)

        dtempij = -np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (
            rmj @ np.outer(np.cross(rij - rkj, rmj), rmj).T + 
            np.outer(np.cross(rij - rkj, rmj), rmj) @ rmj.T
        )
        dij = -np.outer(rmj,rkj)/ (np.dot(rmj.T, rmj) * np.linalg.norm(rkj)) + dtempij
        dtempik = np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (
            rmj @ np.outer(np.cross(rij, rmj), rmj).T + 
            np.outer(np.cross(rij, rmj), rmj) @ rmj.T
        )
        dik = np.outer(rmj,rkj) / (np.dot(rmj.T, rmj) * np.linalg.norm(rkj)) + dtempik
        dil = np.zeros((3, 3))

        dll = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (
            rnk @ np.outer(np.cross(rkj, rnk), rnk).T + 
            np.outer(np.cross(rkj, rnk), rnk) @ rnk.T
        )
        dtemplk = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (
            rnk @ np.outer(np.cross(rkl - rkj, rnk), rnk).T + 
            np.outer(np.cross(rkl - rkj, rnk), rnk) @ rnk.T
        )
        dlk = -rnk @ rkj.T / (np.dot(rnk.T, rnk) * np.linalg.norm(rkj)) + dtemplk
        dtemplj = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (
            rnk @ np.outer(np.cross(rnk, rkl), rnk).T + 
            np.outer(np.cross(rnk, rkl), rnk) @ rnk.T
        )
        dlj = -np.outer(rnk,rkj) / (np.dot(rnk.T, rnk) * np.linalg.norm(rkj)) + dtemplj

        dT1jj = 1 / np.dot(rkj.T, rkj) * ((-1 + 2 * np.dot(rij.T, rkj) / np.dot(rkj.T, rkj)) * rkj - rij)
        dT2jj = 1 / np.dot(rkj.T, rkj) * (2 * np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * rkj - rkl)
        djj = np.outer(di,dT1jj) + (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * dij - \
            (np.outer(di,dT2jj) + np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dlj)

        dT1jk = 1 / np.dot(rkj.T, rkj) * (-2 * np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * rkj + rij)
        dT2jk = 1 / np.dot(rkj.T, rkj) * ((1 - 2 * np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj)) * rkj + rkl)
        djk = np.outer(di,dT1jk) + (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * dik - \
            (np.outer(dl,dT2jk) + np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dlk)


        dT1kk = dT2jk
        dT2kk = dT1jk
        dkk = np.outer(dl,dT1kk) + (np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) - 1) * dlk - (di @ dT2kk.T + np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * dik)

        # Hp = np.block([
        #     [djj, djk, dij.T, dlj.T],
        #     [djk.T, dkk, dik.T, dlk.T],
        #     [dij, dik, dii, dil],
        #     [dlj, dlk, dil.T, dll]
        # ])

        Khe = np.zeros((12, 12))    
        # (Kspr * (Jhe @ Jhe.T) + Rspr * Hp)

    return he, Rhe, Khe

# def FoldKe(Cood, List, kpi, h0=None, L0=None, CM=None):
#     rkj = (Cood[List[1]] - Cood[List[0]]).T
#     rij = (Cood[List[2]] - Cood[List[0]]).T
#     rkl = (Cood[List[1]] - Cood[List[3]]).T
    
#     rmj = np.cross(rij, rkj)
#     rnk = np.cross(rkj, rkl)
    
#     sgn = ((np.abs(np.dot(rnk, rij)) > 1e-8) * np.sign(np.dot(rnk, rij)) + (np.abs(np.dot(rnk, rij)) <= 1e-8) * 1)
#     he = np.real(np.arccos(np.dot(rmj, rnk) / (np.linalg.norm(rmj) * np.linalg.norm(rnk))))
#     he = np.real(sgn * he)
    
#     if he < 0:
#         he = 2 * np.pi + he
    
#     Rspr, Kspr = CM(he, h0, kpi, L0)
        
#     di = np.linalg.norm(rkj) / (np.dot(rmj, rmj)) * rmj
#     dl = -np.linalg.norm(rkj) / (np.dot(rnk, rnk)) * rnk
#     dj = (np.dot(rij, rkj) / (np.dot(rkj, rkj)) - 1) * di - np.dot(rkl, rkj) / (np.dot(rkj, rkj)) * dl
#     dk = -np.dot(rij, rkj) / (np.dot(rkj, rkj)) * di + (np.dot(rkl, rkj) / (np.dot(rkj, rkj)) - 1) * dl
    
#     Jhe = np.array([dj, dk, di, dl])
#     Rhe = Rspr @ Jhe
    
#     dii = -np.linalg.norm(rkj) / (np.dot(rmj, rmj) ** 2) * ((rmj @ np.cross(rkj, rmj).T) + (rmj @ np.cross(rkj, rmj).T).T)
    
#     dtempij = -np.linalg.norm(rkj) / (np.dot(rmj, rmj) ** 2) * (rmj @ (np.cross(rij - rkj, rmj).T) + (np.cross(rij - rkj, rmj)) @ rmj.T)
#     dij = -rmj @ rkj.T / (np.dot(rmj, rmj) * np.linalg.norm(rkj)) + dtempij
    
#     dtempik = np.linalg.norm(rkj) / (np.dot(rmj, rmj) ** 2) * (rmj @ (np.cross(rij, rmj).T) + (np.cross(rij, rmj)) @ rmj.T)
#     dik = rmj @ rkj.T / (np.dot(rmj, rmj) * np.linalg.norm(rkj)) + dtempik
    
#     dil = np.zeros((3, 3))
    
#     dll = np.linalg.norm(rkj) / (np.dot(rnk, rnk) ** 2) * (rnk @ np.cross(rkj, rnk).T + (rnk @ np.cross(rkj, rnk).T).T)
    
#     dtemplk = np.linalg.norm(rkj) / (np.dot(rnk, rnk) ** 2) * (rnk @ (np.cross(rkl - rkj, rnk).T) + (np.cross(rkl - rkj, rnk)) @ rnk.T)
#     dlk = -rnk @ rkj.T / (np.dot(rnk, rnk) * np.linalg.norm(rkj)) + dtemplk
    
#     dtemplj = np.linalg.norm(rkj) / (np.dot(rnk, rnk) ** 2) * (rnk @ (np.cross(rnk, rkl).T) + (rnk @ np.cross(rnk, rkl).T).T)
#     dlj = rnk @ rkj.T / (np.dot(rnk, rnk) * np.linalg.norm(rkj)) + dtemplj
    
#     dT1jj = 1 / (np.dot(rkj, rkj)) * ((-1 + 2 * np.dot(rij, rkj) / (np.dot(rkj, rkj))) * rkj - rij)
#     dT2jj = 1 / (np.dot(rkj, rkj)) * (2 * np.dot(rkl, rkj) / (np.dot(rkj, rkj)) * rkj - rkl)
#     djj = di @ dT1jj.T + (np.dot(rij, rkj) / (np.dot(rkj, rkj)) - 1) * dij - (dl @ dT2jj.T + np.dot(rkl, rkj) / (np.dot(rkj, rkj)) * dlj)
    
#     dT1jk = 1 / (np.dot(rkj, rkj)) * (-2 * np.dot(rij, rkj) / (np.dot(rkj, rkj)) * rkj + rij)
#     dT2jk = 1 / (np.dot(rkj, rkj)) * ((1 - 2 * np.dot(rkl, rkj) / (np.dot(rkj, rkj))) * rkj + rkl)
#     djk = di @ dT1jk.T + (np.dot(rij, rkj) / (np.dot(rkj, rkj)) - 1) * dik - (dl @ dT2jk.T + np.dot(rkl, rkj) / (np.dot(rkj, rkj)) * dlk)
    
#     dT1kk = dT2jk
#     dT2kk = dT1jk
#     dkk = dl @ dT1kk.T + (np.dot(rkl, rkj) / (np.dot(rkj, rkj)) - 1) * dlk - (di @ dT2kk.T + np.dot(rij, rkj) / (np.dot(rkj, rkj)) * dik)
    
#     Hp = np.array([[djj, djk, dij.T, dlj.T],
#                     [djk, dkk, dik.T, dlk.T],
#                     [dij, dik, dii, dil],
#                     [dlj, dlk, dil.T, dll]])
    
#     Khe = Kspr * (Jhe @ Jhe.T) + Rspr @ Hp
    
#     return he, Rhe, Khe

