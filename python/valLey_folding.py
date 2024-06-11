# =========== VALLEY FOLDING =========================================== 

# Define geometry and material 

import numpy as np
from ConfigMiura import ConfigMiura
from Ogden import Ogden
from PrepareData import PrepareData
from PathAnalysis import PathAnalysis
from EnhancedLinear import EnhancedLinear
from VisualFold import VisualFold
from PostProcess import PostProcess

Node = np.array([[0,0,0],
                [ 1, 0, 0],
                [1, 1, 0.1],
                [0, 1, 0]])
Panel = [[0, 1, 3],
    [1, 2, 3]]
BDRY = np.array([[0, 1],
                [0, 3]])

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
    [1, 0, 1, 1],
    [3, 1, 0, 1] ])

indp = [2]
ff = np.ones(len(indp))
Load = np.column_stack((indp, np.zeros_like(indp), np.zeros_like(indp), ff))
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
instdof = -(indp[0] * 3 - 2)
# interv = 0
# endicrm = U_his.shape[1]
VisualFold(U_his, truss, angles, 'none', 'valley_folding', 0.05, LF_his, instdof, [-np.inf, np.inf, -np.inf, np.inf])
 
# If do not need load-displacement diagram:
# VisualFold(U_his[:, :interv:endicrm], truss, angles, 'none', 'miura5x5fold', 0.0001)

# Plot stored energy vs. pseudo time
# Red line is the total profile.  Between red and cyan is the folding
# energy. Between cyan and magenta is the portion of energy for bending. 
# Below magenta is the stretching energy of bars.
# STAT = PostProcess(Data, truss, angles)
# plt.figure()
# plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.PE.strain, 'r-', linewidth=2)  # Total profile
# plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.bend.PE + STAT.bar.PE, 'c-')  # Folding energy
# plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.bar.PE, 'm-')  # Stretching energy of bars
# plt.xlabel('Increment Number (Pseudo-time)', fontsize=14)
# plt.ylabel('Stored Energy', fontsize=14)
# plt.grid(True)
# plt.show()



    # plt.clf()
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax = plot_panels(Node, Panel, ax)
#plt.draw()
#plt.pause(0.01)
