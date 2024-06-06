import numpy as np

def BarKe(u, B, L, CM, A):
    du = u[:3] - u[3:6]
    Du = np.concatenate((du, -du))
    
    Ex = (B @ u / L + 0.5 * (du.T @ du) / L**2)
    Sx, Et, Wb = CM(Ex)
    Fx = Sx * A

    Rbe = np.dot(Fx, (B.A + Du / L))

    Kel = np.dot(B.T, B)
    Kg = Fx / L * np.block([[np.eye(3), -np.eye(3)], [-np.eye(3), np.eye(3)]])
    K1 = ((Du*B.A.T) + (Du*B.A.T).T) / L
    K2 = (np.outer(Du, Du.T)) / L**2
    Kbe = (Et * A / L)[0] * (Kel + K1 + K2) + Kg
    return Ex, Rbe, Kbe
    
    # Rbe = Fx*(B + Du/L)

    # Kel = B.T @ B
    # Kg = Fx / L * np.block([[np.eye(3), -np.eye(3)], [-np.eye(3), np.eye(3)]])
    # K1 = ((np.dot(Du, B) + np.dot(Du, B).T))/L
    # K2 = np.dot(Du, Du.T)/L**2
    # Kbe = np.dot(Et, A)/L * (Kel + K1 + K2) + Kg
    
    # return Ex, Rbe, Kbe

