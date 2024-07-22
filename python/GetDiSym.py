import numpy as np

def GetDiSym(N, h, lyr, phi):
    '''The function `GetDiSym` generates node coordinates and panel connectivity for a symmetric geometry
    based on input parameters.
    
    Parameters
    ----------
    N
        N is the number of nodes in each layer.
    h
        The parameter `h` in the `GetDiSym` function represents the height of the structure. It is used to
    calculate the Z coordinates of the nodes in the structure based on the number of layers (`lyr`) and
    the number of nodes per layer (`N`).
    lyr
        The parameter `lyr` in the `GetDiSym` function represents the number of layers in the structure. It
    is used to determine the number of layers in the structure and to perform calculations based on the
    specified number of layers.
    phi
        The `phi` parameter in the `GetDiSym` function represents the rotation angle for each layer in the
    structure. It can be a single value if the structure has the same rotation angle for all layers, or
    it can be an array of rotation angles if each layer has a different rotation angle.
    
    Returns
    -------
        The function `GetDiSym` returns two arrays: `Node` and `Panel`.
    
    '''
    if len(np.atleast_1d(phi)) == 1:
        rotangle = (np.arange(lyr)) * phi
    else:
        rotangle = phi
    
    rdl = np.zeros((lyr, N))
    for i in range(lyr):
        rdl[i, :] = np.linspace(rotangle[i], 2 * np.pi / N * (N - 1) + rotangle[i], N)
    
    Xcood = np.cos(rdl.T.reshape(-1, 1))
    Ycood = np.sin(rdl.T.reshape(-1, 1))
    Zcood = h * np.repeat(np.arange(lyr), N)
    Node = np.column_stack((Xcood, Ycood, Zcood))
    
    PMat = np.zeros((2 * N * (lyr - 1), 3), dtype=int)
    for i in range(lyr - 1):
        PMat[(i * 2 * N):(i * 2 * N + N), :] = np.column_stack((np.arange(1, N + 1), np.arange(1, N + 1) + N, np.mod(np.arange(1, N + 1), N) + N + 1)) + i * N
        PMat[(i * 2 * N + N):(i * 2 * N + 2 * N), :] = np.column_stack((np.arange(1, N + 1), np.mod(np.arange(1, N + 1), N) + 1, np.mod(np.arange(1, N + 1), N) + N + 1)) + i * N
    
    Panel = [PMat[i, :] for i in range(PMat.shape[0])]
    
    return Node, Panel

