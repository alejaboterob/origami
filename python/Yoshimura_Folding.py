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

from ConfigYoshimura2 import ConfigYoshimura

from PostProcess import PostProcess

from matplotlib import pyplot as plt


# YOSHIMURA FOLDING

# Define geometry and material parameters


# Example usage
sec_hor = 8  # Number of sections around the circumference
sec_vert = 2  # Number of sections along the height
radius = 1  # Radius of the cylinder
height = 2  # Height of the cylinder
fdang = 0.0  # Folding angle (for future use or deformation)

Node, Panel, BDRY = ConfigYoshimura(sec_hor, sec_vert, radius, height, fdang)


# Call ConfigMiura to get Node and Panel arrays
# Node, Panel, BDRY = ConfigMiura(alpha, n, l, m)

# Node, Panel, BDRY = ConfigMiura(sec_hor, sec_vert, theta, a, b, fdang)
# Node, Panel, _ = ConfigMiura(sec_hor, sec_vert, theta, a, b, fdang)

# Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
# E0-stretching stiffness; Abar-bar area (assume uniform)
Kf = 1e-1
Kb = Kf * 1e5
E0 = 1e6
Abar = 1e-1

# Left and right limits for the linear range of rotational stiffness
limlft = 0.1
limrht = 360 - 0.1

# Define material constitutive and rotational spring functions
BarMater = lambda Ex: Ogden(Ex, E0)  # Define bar material constitutive
RotSpring = lambda he, h0, kpi, L0: EnhancedLinear(he, h0, kpi, L0, limlft, limrht)

# Maximum increment number & initial load factor
MaxIcr = 60
blam = 0.5

# Node indices on the left side (along x-axis)
leftx = np.arange(0, (2 * sec_vert + 1))

# Node indices on the left side (along z-axis)
leftz = np.arange(0, (2 * sec_vert + 1) + 1, 2)

# Node indices on the right side (along z-axis)
rightz = np.arange(0, (2 * sec_vert + 1) + 1, 2) + (2 * sec_vert + 1) * (2 * sec_hor)

# Node indices on the right side (along x-axis)
rightxp = np.arange(1, (2 * sec_vert + 1), 2) + (2 * sec_vert + 1) * (2 * sec_hor)

# Set up the support conditions
Supp = np.array([[0, 0, 1, 0],  # Example of a fixed node (node 0, z direction constrained)
                    *zip(leftx, np.ones_like(leftx), np.zeros_like(leftx), np.zeros_like(leftx)),
                    *zip(leftz, np.zeros_like(leftz), np.zeros_like(leftz), np.ones_like(leftz)),
                    *zip(rightz, np.zeros_like(rightz), np.zeros_like(rightz), np.ones_like(rightz))])


# indp = np.arange(0, (sec_vert * 2 + 1)) + (sec_vert * 2 + 1) * (sec_hor * 2)   # Revisar -1 
# ff = -np.ones(len(indp))
# Load = np.column_stack((indp, ff, np.zeros_like(indp), np.zeros_like(indp)))
# indp = Load[:, 0]

Load = np.array([[0, 0, 0, 1],
                [1, 0, 0, 1],
                [2, 0, 0, 1],
                [3, 0, 0, 1]])

# Perform analysis
truss, angles, F = PrepareData(Node, Panel, Supp, Load, BarMater, RotSpring, Kf, Kb, Abar)
truss['U0'] = np.zeros(3 * truss['Node'].shape[0])

U_his, LF_his, Data = PathAnalysis(truss, angles, F, blam, MaxIcr)
U_his = np.real(U_his)
LF_his = np.real(LF_his)

print(truss)
print(angles)
print(F)

indp = Load[:, 0]


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



