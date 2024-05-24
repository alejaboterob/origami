import numpy as np

def icross(a, b):
    return np.cross(a, b, axis=0)

def fold_ke(Cood, List, kpi, h0=None, L0=None, CM=None):
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

        if nargout > 1:
            di = np.linalg.norm(rkj) / np.dot(rmj.T, rmj) * rmj
            dl = -np.linalg.norm(rkj) / np.dot(rnk.T, rnk) * rnk
            dj = (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * di - np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dl
            dk = -np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * di + (np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) - 1) * dl

            Jhe = np.vstack([dj, dk, di, dl])
            Rhe = Rspr * Jhe

        if nargout > 2:
            dii = -np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * ((rmj @ icross(rkj, rmj).T) + (rmj @ icross(rkj, rmj).T).T)
            dtempij = -np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (rmj @ (icross(rij - rkj, rmj)).T + icross(rij - rkj, rmj) @ rmj.T)
            dij = -rmj @ rkj.T / (np.dot(rmj.T, rmj) * np.linalg.norm(rkj)) + dtempij
            dtempik = np.linalg.norm(rkj) / (np.dot(rmj.T, rmj)**2) * (rmj @ (icross(rij, rmj)).T + icross(rij, rmj) @ rmj.T)
            dik = rmj @ rkj.T / (np.dot(rmj.T, rmj) * np.linalg.norm(rkj)) + dtempik
            dil = np.zeros((3, 3))

            dll = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (rnk @ icross(rkj, rnk).T + (rnk @ icross(rkj, rnk).T).T)
            dtemplk = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (rnk @ (icross(rkl - rkj, rnk)).T + icross(rkl - rkj, rnk) @ rnk.T)
            dlk = -rnk @ rkj.T / (np.dot(rnk.T, rnk) * np.linalg.norm(rkj)) + dtemplk
            dtemplj = np.linalg.norm(rkj) / (np.dot(rnk.T, rnk)**2) * (rnk @ (icross(rnk, rkl)).T + (rnk @ icross(rnk, rkl).T).T)
            dlj = rnk @ rkj.T / (np.dot(rnk.T, rnk) * np.linalg.norm(rkj)) + dtemplj

            dT1jj = 1 / (np.dot(rkj.T, rkj)) * ((-1 + 2 * np.dot(rij.T, rkj) / np.dot(rkj.T, rkj)) * rkj - rij)
            dT2jj = 1 / (np.dot(rkj.T, rkj)) * (2 * np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * rkj - rkl)
            djj = di @ dT1jj.T + (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * dij - (dl @ dT2jj.T + np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dlj)

            dT1jk = 1 / (np.dot(rkj.T, rkj)) * (-2 * np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * rkj + rij)
            dT2jk = 1 / (np.dot(rkj.T, rkj)) * ((1 - 2 * np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj)) * rkj + rkl)
            djk = di @ dT1jk.T + (np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) - 1) * dik - (dl @ dT2jk.T + np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) * dlk)

            dT1kk = dT2jk
            dT2kk = dT1jk
            dkk = dl @ dT1kk.T + (np.dot(rkl.T, rkj) / np.dot(rkj.T, rkj) - 1) * dlk - (di @ dT2kk.T + np.dot(rij.T, rkj) / np.dot(rkj.T, rkj) * dik)

            Hp = np.block([
                [djj, djk, dij.T, dlj.T],
                [djk.T, dkk, dik.T, dlk.T],
                [dij, dik, dii, dil],
                [dlj, dlk, dil.T, dll]
            ])

            Khe = (Kspr * (Jhe @ Jhe.T) + Rspr * Hp)

    return he, Rhe, Khe

# Placeholder CM function (to be implemented according to the specific model)
def CM(he, h0, kpi, L0):
    # Implement the constitutive model function here
    Rspr = ...  # Example placeholder
    Kspr = ...  # Example placeholder
    return Rspr, Kspr

# Placeholder nargout function for compatibility
def nargout():
    import inspect
    frame = inspect.currentframe().f_back
    return frame.f_globals['nargout'] if 'nargout' in frame.f_globals else 1
