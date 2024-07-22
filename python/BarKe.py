import numpy as np


def BarKe(u, B, L, CM, A):
    '''The function `BarKe` calculates the strain, reaction forces, and stiffness matrix for a bar element
    in a structure.

    Parameters
    ----------
    u
        It seems like the description of the parameters is missing. Could you please provide more
    information about the parameters `B`, `L`, `CM`, and `A` as well? This will help me provide a more
    accurate explanation of the function `BarKe`.
    B
        The parameter B is typically a matrix that represents the shape functions of an element in a finite
    element analysis. It is used to interpolate the field variables within the element based on the
    nodal values.
    L
        The parameter L in the function `BarKe` represents the length of the bar element. It is used in
    various calculations within the function to determine the strain, stress, forces, and stiffness of
    the bar element.
    CM
        It seems like the definition of the `CM` function is missing in the provided code snippet. The `CM`
    function is used in the `BarKe` function, but its implementation is not shown.
    A
        The parameter A in the function `BarKe` represents the cross-sectional area of the bar element. It
    is used in the calculation of the force Fx and the stiffness matrix Kbe.

    Returns
    -------
        The function `BarKe` returns three values: `Ex`, `Rbe`, and `Kbe`. `Ex` represents the strain in
    the bar element, `Rbe` represents the internal forces in the bar element, and `Kbe` represents the
    stiffness matrix of the bar element.

    '''
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

