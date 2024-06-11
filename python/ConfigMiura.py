import numpy as np

def ConfigMiura(x_divs, y_divs, theta, a, b, gmma):
    
    theta = np.pi * theta / 180
    gmma = gmma * np.pi / 180
    numx = 2 * x_divs
    numy = 2 * y_divs
    h = a * np.sin(theta) * np.sin(gmma)
    s = b * np.cos(gmma) * np.tan(theta) / np.sqrt(1 + np.cos(gmma) ** 2 * np.tan(theta) ** 2)
    l = a * np.sqrt(1 - np.sin(gmma) ** 2 * np.sin(theta) ** 2)
    v = b / np.sqrt(1 + np.cos(gmma) ** 2 * np.tan(theta) ** 2)
    X = np.repeat(s * np.arange(numx + 1), numy + 1).reshape(numy + 1, numx + 1)
    Y = np.repeat(l * np.arange(numy + 1), numx + 1).reshape(numy + 1, numx + 1)
    Y[:, 1::2] = Y[:, 1::2] + v
    Z = np.zeros_like(X)
    Z[1::2, :] = h

    # Reshape each array to a column vector
    X_reshaped = X.reshape(-1, 1)  # WHYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
    Y_reshaped = Y.T.reshape(-1, 1)
    Z_reshaped = Z.T.reshape(-1, 1)

    # Concatenate the reshaped arrays column-wise
    NODE = np.hstack((X_reshaped, Y_reshaped, Z_reshaped))

    BDRY = np.zeros((2 * numx + 2 * numy, 2), dtype=int)
    Node_idx = np.arange(1, len(NODE) + 1).reshape(numy + 1, numx + 1)
    Lbdry = Node_idx[0, :]
    Rbdry = Node_idx[-1, :]
    Bbdry = Node_idx[:, 0]
    Tbdry = Node_idx[:, -1]

    # Specify boundary
    BDRY[:len(Lbdry) - 1, :] = np.column_stack((Lbdry[:-1], Lbdry[1:]))
    count = len(Lbdry) - 1
    BDRY[count:count + len(Bbdry) - 1, :] = np.column_stack((Bbdry[:-1], Bbdry[1:]))
    count += len(Bbdry) - 1
    BDRY[count:count + len(Rbdry) - 1, :] = np.column_stack((Rbdry[:-1], Rbdry[1:]))
    count += len(Rbdry) - 1
    BDRY[count:count + len(Tbdry) - 1, :] = np.column_stack((Tbdry[:-1], Tbdry[1:]))

    BDRY = BDRY - 1 # OJO 

    # PANEL = [None] * (numx * numy)
    PANEL = np.zeros((numx * numy, 4), dtype=int)
    k = 0
    for j in range(numy):
        for i in range(numx):
            n1 = (i) * (numy + 1) + j   # OJO +1 igual que matlab
            n2 = (i+1) * (numy + 1) + j   # OJO 
            PANEL[k] = [n1, n2, n2 + 1, n1 + 1]
            k += 1

    return NODE, PANEL, BDRY




# import numpy as np

# def ConfigMiura(x_divs, y_divs, theta, a, b, gmma):
#     theta = np.pi * theta / 180
#     gmma = gmma * np.pi / 180
#     numx = 2 * x_divs
#     numy = 2 * y_divs
#     h = a * np.sin(theta) * np.sin(gmma)
#     s = b * np.cos(gmma) * np.tan(theta) / np.sqrt(1 + np.cos(gmma) ** 2 * np.tan(theta) ** 2)
#     l = a * np.sqrt(1 - np.sin(gmma) ** 2 * np.sin(theta) ** 2)
#     v = b / np.sqrt(1 + np.cos(gmma) ** 2 * np.tan(theta) ** 2)
#     X = np.tile(s * np.arange(0, numx + 1), (numy + 1, 1))
#     Y = np.tile(l * np.arange(0, numy + 1)[:, np.newaxis], (1, numx + 1))
#     Y[:, 1::2] += v
#     Z = np.zeros_like(X)
#     Z[1::2, :] = h
#     NODE = np.column_stack((X.flatten(), Y.flatten(), Z.flatten()))
    
#     # Specify boundary

#     BDRY = np.zeros((2 * numx + 2 * numy, 2))
#     Node_idx = np.reshape(np.arange(0, NODE.shape[0]), X.shape)
#     Lbdry = Node_idx[:, 0]
#     Rbdry = Node_idx[:, -1]
#     Bbdry = Node_idx[0, :]
#     Tbdry = Node_idx[-1, :]

#     BDRY[0:len(Lbdry) - 2, :] = np.column_stack((Lbdry[0:-2], Lbdry[1:]))
#     count = len(Lbdry) - 1
#     BDRY[count:count + len(Bbdry) - 2, :] = np.column_stack((Bbdry[:-2], Bbdry[1:]))
#     count += len(Bbdry) - 1
#     BDRY[count:count + len(Rbdry) - 2, :] = np.column_stack((Rbdry[:-2], Rbdry[1:]))
#     count += len(Rbdry) - 1
#     BDRY[count:count + len(Tbdry) - 2, :] = np.column_stack((Tbdry[:-2], Tbdry[1:]))
#     count += len(Tbdry) - 1

#     k = 0
#     PANEL = []
#     for j in range(0, numy):
#         for i in range(0, numx):
#             k += 1
#             n1 = i * (numy + 2) + j
#             n2 = (i-1) * (numy + 2) + j
#             PANEL.append([n1, n2, n2 + 1, n1 + 1])

#     return NODE, PANEL, BDRY
