import numpy as np

def ConfigYoshimura(alpha, n, l, m):

    '''The function `ConfigMiura` generates node coordinates, panel connectivity, and boundary information
    for a Miura pattern based on input parameters.

    Parameters
    ----------
    x_divs
        X divisions in the Miura pattern
    y_divs
        The `y_divs` parameter in the `ConfigMiura` function represents the number of divisions along the
    y-axis in the Miura pattern configuration. It determines how many segments the pattern will have
    vertically. Increasing `y_divs` will result in more detailed and smaller segments along the y-axis
    theta
        Theta is the angle in degrees that is converted to radians for calculations in the function.
    a
        Parameter 'a' represents the length of the Miura fold along the x-axis.
    b
        Parameter `b` represents the length of the Miura fold unit cell in the y-direction. It is used in
    the Miura configuration calculation to determine the geometry of the fold pattern.
    gmma
        The parameter `gmma` in the `ConfigMiura` function represents the angle `gamma` in degrees. It is
    converted to radians within the function using the formula `gmma = gmma * np.pi / 180`. This
    conversion allows the angle to be used in trigonometric calculations

    Returns
    -------
        The function `ConfigMiura` returns three arrays: `NODE`, `PANEL`, and `BDRY`.

    '''

    # Setting some parameters
    openingAngle =30 *np.pi/180 # opening angle in degrees
    halfLayerheight =1 # height of half a layer
    eps =0.5 # length of the small edges
    # length of " diagonal " creases
    diagLen = halfLayerheight /  np.sin( openingAngle )
    # length of longer horizontal creases
    horzLen =2* halfLayerheight / np.tan( openingAngle ) + eps
    # size of the pattern
    numDoubleNodes =2
    numLayers =4
    









    beta = (180-((4*n-2)*180)/4*n)/2 * np.pi/180 # Introduction angle
    theta = 2*np.arccos((1/np.tan(alpha))*np.tan(beta)) # Dihedral angle between two adjacent faces during folding 
    H = l / (4*np.cos(beta)*np.sin(np.pi/(4*n)))  # Height of the entire structure
    S = 2*H  # Span of the entire structure
    L = m * (l/2) * np.tan(alpha)*np.sin(theta/2)  # Length of the entire structure





    # 2 * np.arccos(1/((1/np.tan(beta))*np.tan(alpha)))  

    # Node coordinates list
    nodes = []
    
    # Calculate coordinates for each node
    for i in range(m + 1):
        for j in range(n + 1):
            x = i * (l / 2) * np.tan(alpha) * np.sin(theta / 2)
            y = j * l / (4 * np.cos(beta))
            z = (-1)**(i + j) * H / 2  # Alternate z values to create the zigzag pattern
            nodes.append((x, y, z))


    nodes = []

    # Define base points
    x0, y0, z0 = 0, 0, 0
    x1, y1, z1 = l * np.cos(alpha), 0, l * np.sin(alpha)
    x2, y2, z2 = l * np.cos(alpha) / 2, l * np.sin(alpha) / 2, H
    x3, y3, z3 = -l * np.cos(alpha) / 2, l * np.sin(alpha) / 2, H
    x4, y4, z4 = -l * np.cos(alpha), 0, l * np.sin(alpha)
    x5, y5, z5 = -l * np.cos(alpha) / 2, -l * np.sin(alpha) / 2, H
    x6, y6, z6 = l * np.cos(alpha) / 2, -l * np.sin(alpha) / 2, H

    # Append the nodes to the list
    nodes.append((x0, y0, z0))
    nodes.append((x1, y1, z1))
    nodes.append((x2, y2, z2))
    nodes.append((x3, y3, z3))
    nodes.append((x4, y4, z4))
    nodes.append((x5, y5, z5))
    nodes.append((x6, y6, z6))






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




import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def yoshimura_pattern(n_units, n_layers, L, beta_deg):
    nodes = []
    beta_rad = np.radians(beta_deg)
    
    for i in range(n_layers + 1):
        for j in range(n_units + 1):
            x = j * L
            y = i * L * np.cos(beta_rad)
            z = (i % 2) * L * np.sin(beta_rad) * (1 if j % 2 == 0 else -1)
            nodes.append((x, y, z))
    
    return nodes

def yoshimura_panels(n_units, n_layers):
    panels = []
    
    for i in range(n_layers):
        for j in range(n_units):
            current = i * (n_units + 1) + j
            next_row = (i + 1) * (n_units + 1) + j
            
            if i % 2 == 0:
                if j % 2 == 0:
                    panels.append([current, current + 1, next_row + 1])
                    panels.append([current, next_row + 1, next_row])
                else:
                    panels.append([current, current + 1, next_row])
                    panels.append([current + 1, next_row + 1, next_row])
            else:
                if j % 2 == 0:
                    panels.append([current, current + 1, next_row])
                    panels.append([current + 1, next_row + 1, next_row])
                else:
                    panels.append([current, current + 1, next_row + 1])
                    panels.append([current, next_row + 1, next_row])
    
    return panels

def plot_yoshimura(nodes, panels):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xs, ys, zs = zip(*nodes)
    ax.scatter(xs, ys, zs, c='r', marker='o')

    for i, (x, y, z) in enumerate(nodes):
        ax.text(x, y, z, f'{i+1}', color='blue')
    
    for panel in panels:
        verts = [nodes[panel[0]], nodes[panel[1]], nodes[panel[2]]]
        tri = Poly3DCollection([verts], color='cyan', edgecolor='k')
        ax.add_collection3d(tri)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()

def fold_pattern(nodes, fold_amount):
    folded_nodes = []
    for (x, y, z) in nodes:
        z_folded = z * fold_amount
        folded_nodes.append((x, y, z_folded))
    return folded_nodes

# Example usage
n_units = 6      # Number of units horizontally
n_layers = 4     # Number of layers vertically
L = 100          # Length of each unit
beta_deg = 30    # Beta angle in degrees
fold_amount = 0.5  # Folding amount factor (0: flat, 1: fully folded)

nodes = yoshimura_pattern(n_units, n_layers, L, beta_deg)
panels = yoshimura_panels(n_units, n_layers)
folded_nodes = fold_pattern(nodes, fold_amount)

print("Node coordinates for Yoshimura pattern:")
for i, (x, y, z) in enumerate(nodes):
    print(f"Node {i+1}: x = {x:.2f}, y = {y:.2f}, z = {z:.2f}")

plot_yoshimura(folded_nodes, panels)
