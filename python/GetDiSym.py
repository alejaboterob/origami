import numpy as np

def GetDiSym(N, h, lyr, phi):
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

