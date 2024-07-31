# This Python script is a simulation and analysis program for Miura folding, a type of origami folding
# pattern. Here is a breakdown of what the script does:


import numpy as np
from ConfigMiura import ConfigMiura
from Ogden import Ogden
from PrepareData import PrepareData
from PathAnalysis import PathAnalysis
from EnhancedLinear import EnhancedLinear
from VisualFold import VisualFold
from displacement import displacement
from GraphPostProcess import GraphPostProcess
from PlotFinalConfig import PlotFinalConfig

from PostProcess import PostProcess

from matplotlib import pyplot as plt


# YOSHIMURA FOLDING

# Define geometry and material parameters

alpha = np.pi/4 # Plane angle between the valley and mountain folds in the crease scheme  The plane angle α must be 0 < α ≤ π/4 for the pattern to be folded
l = 1  # Length of the long side of the isosceles triangle
n = 1  # Number of complete valley folds in the crease scheme
m = 1  # Number of repetitions of the basic crease scheme



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def yoshimura_pattern(n_units, n_layers, L, beta):
    nodes = []
    
    for i in range(n_layers + 1):
        for j in range(n_units + 1):
            x = j * L
            y = i * L * np.cos(beta)
            z = (i % 2) * L * np.sin(beta) * (1 if j % 2 == 0 else -1)
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
beta = 30  *np.pi/180  # Beta angle in degrees
fold_amount = 0.5  # Folding amount factor (0: flat, 1: fully folded)

Node = np.array(yoshimura_pattern(n_units, n_layers, L, beta))
Panel = np.array(yoshimura_panels(n_units, n_layers))
folded_nodes = fold_pattern(Node, fold_amount)

print("Node coordinates for Yoshimura pattern:")
for i, (x, y, z) in enumerate(Node):
    print(f"Node {i+1}: x = {x:.2f}, y = {y:.2f}, z = {z:.2f}")

plot_yoshimura(folded_nodes, Panel)






sec_hor = 5  # Number of unit cells in horizontal direction
sec_vert = 5  # Number of unit cells in vertical direction
# Geometry of the Miura: a, b are the edge lengths of each parallelogram
# panel; fdang controls the folding angle; theta is the panel angle
theta = 60
a = 2
b = 2
fdang = 15

# Maximum increment number & initial load factor
MaxIcr = 60
blam = 0.5

# Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
# E0-stretching stiffness; Abar-bar area (assume uniform)
Kf = 1e-1
Kb = Kf * 1e5
E0 = 1e6
Abar = 1e-1

# Left and right limits for the linear range of rotational stiffness
limlft = 0.1
limrht = 360 - 0.1

# Call ConfigMiura to get Node and Panel arrays
# Node, Panel, BDRY = ConfigMiura(alpha, n, l, m)

# Node, Panel, BDRY = ConfigMiura(sec_hor, sec_vert, theta, a, b, fdang)
# Node, Panel, _ = ConfigMiura(sec_hor, sec_vert, theta, a, b, fdang)

# Define material constitutive and rotational spring functions
BarMater = lambda Ex: Ogden(Ex, E0)  # Define bar material constitutive
RotSpring = lambda he, h0, kpi, L0: EnhancedLinear(he, h0, kpi, L0, limlft, limrht)

# Set up boundary conditions
leftx = np.arange(0, (2 * sec_vert + 1))
leftz = np.arange(0, (2 * sec_vert + 1) + 1, 2)
rightz = np.arange(0, (2 * sec_vert + 1) + 1, 2) + (2 * sec_vert + 1) * (2 * sec_hor)
rightxp = np.arange(1, (2 * sec_vert + 1), 2) + (2 * sec_vert + 1) * (2 * sec_hor)

Supp = np.array([[0, 0, 1, 0],
                 *zip(leftx, np.ones_like(leftx), np.zeros_like(leftx), np.zeros_like(leftx)),
                 *zip(leftz, np.zeros_like(leftz), np.zeros_like(leftz), np.ones_like(leftz)),
                 *zip(rightz, np.zeros_like(rightz), np.zeros_like(rightz), np.ones_like(rightz))])

indp = np.arange(0, (sec_vert * 2 + 1)) + (sec_vert * 2 + 1) * (sec_hor * 2)   # Revisar -1 
ff = -np.ones(len(indp))
Load = np.column_stack((indp, ff, np.zeros_like(indp), np.zeros_like(indp)))
indp = Load[:, 0]

# Perform analysis
truss, angles, F = PrepareData(Node, Panel, Supp, Load, BarMater, RotSpring, Kf, Kb, Abar)
truss['U0'] = np.zeros(3 * truss['Node'].shape[0])

U_his, LF_his, Data = PathAnalysis(truss, angles, F, blam, MaxIcr)
U_his = np.real(U_his)
LF_his = np.real(LF_his)

print(truss)
print(angles)
print(F)

# Visualize simulation
instdof = -(indp[0] * 3 - 2)
interv = 1
endicrm = U_his.shape[1]
# VisualFold(U_his, truss, angles, 'none', 'valley_folding', 0.05, LF_his, instdof, [-np.inf, np.inf, -np.inf, np.inf])
VisualFold(U_his, truss, angles, LF_his)

# Visualize displacement
instdof = -(indp[0]*3)
displacement(U_his, truss, angles, instdof, 0.1, LF_his)

 
# Plot stored energy vs. pseudo time
# Red line is the total profile.  Between red and cyan is the folding
# energy. Between cyan and magenta is the portion of energy for bending. 
# Below magenta is the stretching energy of bars.
STAT = PostProcess(Data, truss, angles)

GraphPostProcess(U_his, STAT)

PlotFinalConfig(U_his, truss, angles, LF_his)


plt.show()


# print(Node)
# print(Panel)

# from scipy.io import savemat
# variables = {
#     'NODE': Node,
#     'PANEL': Panel,
#     'BDRY': BDRY
# }
# # Save the variables to a .mat file
# savemat('variables.mat', variables)



