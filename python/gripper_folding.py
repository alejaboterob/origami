# =========== GRIPPER FOLDING =========================================== 

# Define geometry and material 

import numpy as np
from ConfigMiura import ConfigMiura
from Ogden import Ogden
from PrepareData import PrepareData
from PathAnalysis import PathAnalysis
from EnhancedLinear import EnhancedLinear
from VisualFold import VisualFold
from PostProcess import PostProcess
from displacement import displacement
from GraphPostProcess import GraphPostProcess
from PlotFinalConfig import PlotFinalConfig
from matplotlib import pyplot as plt

a = 4
b = 2
c = 6
d = 2
e = 2

Node = np.array([[0, -b/2, 0],
                [a, -b/2, 0],
                [0, b/2, 0],
                [a, b/2, -0],

                [(a-e)/2, -(d+b)/2, -d/2],
                [(a+e)/2, -(d+b)/2, -d/2],
                [(a-e)/2, (d+b)/2, -d/2],
                [(a+e)/2, (d+b)/2, -d/2],

                [0, -b/2, c],
                [0, b/2, c],
                [a, -b/2, c],
                [a, b/2, c]])

Panel = [[0, 1, 3, 2],
        [0, 1, 5, 4],
        [2, 3, 7, 6],
        [0, 4, 8],
        [0, 2, 9, 8],
        [2, 6, 9],
        [1, 10, 5],
        [1, 3, 11, 10],
        [3, 7, 11]]

BDRY = np.array([[8, 4],
                [4, 5],
                [5, 10],
                [10, 11],
                [11, 7],
                [7, 6],
                [6, 9],
                [9, 8]])

MaxIcr = 60
blam = 0.01
# Maximum increment number & initial load factor
Kf = 1e-1
Kb = Kf*1e5
E0 = 1e6
Abar = 1e-1
# Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
# E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 0.1
limrht = 180-0.1
# Left and right limits for the linear range of rotational stiffness

# Define material constitutive and rotational spring functions
BarMater = lambda Ex: Ogden(Ex, E0)  # Define bar material constitutive
RotSpring = lambda he, h0, kpi, L0: EnhancedLinear(he, h0, kpi, L0, limlft, limrht)


# Set up boundary conditions
Supp = np.array([[0, 1, 1, 1],
    [1, 1, 1, 1],
    [2, 1, 1, 1],
    [3, 1, 1, 1] ])

indp = [4,5,6,7]
ff = [1,1,-1,-1]* np.ones(len(indp))*3   # OJOOOOOOO  
Load = np.column_stack((indp, np.zeros_like(indp), ff, np.zeros_like(indp)))
indp = Load[:, 0]

# Perform analysis
truss, angles, F = PrepareData(Node, Panel, Supp, Load, BarMater, RotSpring, Kf, Kb, Abar)
truss['U0'] = np.zeros(3 * truss['Node'].shape[0])

U_his, LF_his, Data = PathAnalysis(truss, angles, F, blam, MaxIcr)
U_his = np.real(U_his)
LF_his = np.real(LF_his) 



# from scipy.io import savemat
# variables = {
#     'U_his2': U_his,
#     'LF_his2': LF_his,
#     'Data2': Data

# }
# # Save the variables to a .mat file
# savemat('variables.mat', variables)


# Visualize simulation
instdof = -(indp[0] * 3 +1)
# interv = 0
# endicrm = U_his.shape[1]

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