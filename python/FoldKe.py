import numpy as np

def FoldKe(Cood, List, kpi, h0=None, L0=None, CM=None):
    '''The function `FoldKe` calculates certain geometric properties and returns relevant values based on
    the input parameters.
    
    Parameters
    ----------
    Cood
        It seems like the description of the parameters for the function `FoldKe` is incomplete. Could you
    please provide the descriptions for the following parameters as well?
    List
        The `List` parameter in the `FoldKe` function represents a list of indices that are used to access
    coordinates in the `Cood` array. The function calculates various vectors and angles based on these
    indices.
    kpi
        The `kpi` parameter in the `FoldKe` function seems to be a key parameter used in the calculations
    within the function. It is passed as an argument to the function and is used in various calculations
    involving vectors and matrices. The specific role and meaning of `kpi` would depend on
    h0
        The `h0` parameter in the `FoldKe` function represents the reference angle for the folding
    calculation. It is used in the calculation of the folding energy and related parameters based on the
    input coordinates and list of atoms.
    L0
        In the provided function `FoldKe`, the parameter `L0` is not used directly in the function itself.
    It is passed as an argument but not utilized within the function body.
    CM
        It seems like the `CM` parameter in the `FoldKe` function is a function that takes four arguments
    (`he`, `h0`, `kpi`, `L0`) and returns three values (`Rspr`, `Kspr`, `_`).
    
    Returns
    -------
        The function `FoldKe` returns three values: `he`, `Rhe`, and `Khe`.
    
    '''
    rkj = (Cood[List[1]] - Cood[List[0]]).T
    rij = (Cood[List[2]] - Cood[List[0]]).T
    rkl = (Cood[List[1]] - Cood[List[3]]).T
    rmj = np.cross(rij, rkj)
    rnk = np.cross(rkj, rkl)
    sgn = ((abs(np.dot(rnk.T, rij)) > 1e-8) * np.sign(np.dot(rnk.T, rij)) + 
           (abs(np.dot(rnk.T, rij)) <= 1e-8) * 1)
    he = np.real(np.arccos(np.dot(rmj.T, rnk) / (np.linalg.norm(rmj) * np.linalg.norm(rnk))))
    he = np.real(sgn * he)
    if he < 0:
        he = 2 * np.pi + he

    Rhe = None
    Khe = None  

    if h0 is not None and L0 is not None and CM is not None:
        Rspr, Kspr, _ = CM(he, h0, kpi, L0)
    
        di = np.linalg.norm(rkj) / np.dot(rmj.T, rmj) * rmj
        dl = -np.linalg.norm(rkj) / np.dot(rnk.T, rnk) * rnk
        dj = (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * di - np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dl
        dk = -np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * di + (np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) - 1) * dl

        Jhe = np.vstack([dj, dk, di, dl])
        Rhe = (Rspr * Jhe).flatten()

        dii = -np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (np.outer(rmj, np.cross(rkj, rmj).T) + np.outer(rmj, np.cross(rkj, rmj).T).T)

        dtempij = -np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (
            np.outer(np.cross(rij - rkj, rmj), rmj).T + 
            np.outer(np.cross(rij - rkj, rmj), rmj)
        )
        dij = -np.outer(rmj,rkj)/ (np.dot(rmj.T, rmj) * np.linalg.norm(rkj)) + dtempij
        dtempik = np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (
            np.outer(np.cross(rij, rmj), rmj).T + 
            np.outer(np.cross(rij, rmj), rmj)
        )
        dik = np.outer(rmj,rkj) / (np.dot(rmj.T, rmj) * np.linalg.norm(rkj)) + dtempik
        dil = np.zeros((3, 3))

        dll = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (
            np.outer(np.cross(rkj, rnk), rnk).T + 
            np.outer(np.cross(rkj, rnk), rnk)
        )
        dtemplk = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (
            np.outer(np.cross(rkl - rkj, rnk), rnk).T + 
            np.outer(np.cross(rkl - rkj, rnk), rnk)
        )
        dlk = -np.outer(rnk,rkj) / (np.dot(rnk.T, rnk) * np.linalg.norm(rkj)) + dtemplk
        dtemplj = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (
            np.outer(np.cross(rnk, rkl), rnk).T + 
            np.outer(np.cross(rnk, rkl), rnk)
        )
        dlj = np.outer(rnk,rkj) / (np.dot(rnk.T, rnk) * np.linalg.norm(rkj)) + dtemplj

        dT1jj = 1 / np.dot(rkj.T, rkj) * ((-1 + 2 * np.dot(rij.T, rkj) / np.dot(rkj.T, rkj)) * rkj - rij)
        dT2jj = 1 / np.dot(rkj.T, rkj) * (2 * np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * rkj - rkl)
        djj = np.outer(di,dT1jj) + (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * dij - \
            (np.outer(dl,dT2jj) + np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dlj)

        dT1jk = 1 / np.dot(rkj.T, rkj) * (-2 * np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * rkj + rij)
        dT2jk = 1 / np.dot(rkj.T, rkj) * ((1 - 2 * np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj)) * rkj + rkl)
        djk = np.outer(di,dT1jk) + (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * dik - \
            (np.outer(dl,dT2jk) + np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dlk)


        dT1kk = dT2jk
        dT2kk = dT1jk
        dkk = np.outer(dl,dT1kk) + (np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) - 1) * dlk - (np.outer(di,dT2kk)+ np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * dik)

        Hp = np.block([
            [djj, djk, dij.T, dlj.T],
            [djk.T, dkk, dik.T, dlk.T],
            [dij, dik, dii, dil],
            [dlj, dlk, dil.T, dll]
        ])

        Khe = (Kspr * (np.outer(Jhe.flatten(),Jhe.flatten())) + Rspr * Hp)

        Khe[np.isnan(Khe)] = 0

    return he, Rhe, Khe