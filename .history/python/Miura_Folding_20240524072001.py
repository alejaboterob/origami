import numpy as np
from ConfigMiura import ConfigMiura
from Ogden import Ogden
from PrepareData import PrepareData
from PathAnalysis import PathAnalysis
from EnhancedLinear import EnhancedLinear
from VisualFold import VisualFold
from PostProcess import PostProcess

# MERLIN - Ke Liu, Glaucio H. Paulino
# Ref: K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid origami - An efficient computational approach.'
#      Proceedings of the Royal Society A.
#      K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to capture highly nonlinear behavior of non-rigid origami.'
#      Proceedings of IASS Annual Symposium 2016.

# MIURA FOLDING

# Define geometry and material parameters
sec_hor = 1  # Number of unit cells in horizontal direction
sec_vert = 1  # Number of unit cells in vertical direction
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
Node, Panel, BDRY = ConfigMiura(sec_hor, sec_vert, theta, a, b, fdang)
# Node, Panel, _ = ConfigMiura(sec_hor, sec_vert, theta, a, b, fdang)

# Define material constitutive and rotational spring functions
BarMater = lambda Ex: Ogden(Ex, E0)  # Define bar material constitutive
RotSpring = lambda he, h0, kpi, L0: EnhancedLinear(he, h0, kpi, L0, limlft, limrht)

# Set up boundary conditions
leftx = np.arange(1, (2 * sec_vert + 1) + 1)
leftz = np.arange(1, (2 * sec_vert + 1) + 1, 2)
rightz = np.arange(1, (2 * sec_vert + 1) + 1, 2) + (2 * sec_vert + 1) * (2 * sec_hor)
rightxp = np.arange(2, (2 * sec_vert + 1) + 1, 2) + (2 * sec_vert + 1) * (2 * sec_hor)

Supp = np.array([[1, 0, 1, 0],
                 *zip(leftx, np.ones_like(leftx), np.zeros_like(leftx), np.zeros_like(leftx)),
                 *zip(leftz, np.zeros_like(leftz), np.zeros_like(leftz), np.ones_like(leftz)),
                 *zip(rightz, np.zeros_like(rightz), np.zeros_like(rightz), np.ones_like(rightz))])

indp = np.arange(0, (sec_vert * 2 + 1)) + (sec_vert * 2 + 1) * (sec_hor * 2)+1
ff = -np.ones(len(indp))
Load = np.column_stack((indp-1, ff, np.zeros_like(indp), np.zeros_like(indp)))
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

# # Visualize simulation
# instdof = -(indp[0] * 3 - 2)
# interv = 1
# endicrm = U_his.shape[1]
# VisualFold(U_his[:, :interv:endicrm], truss, angles, 'none', 'miura5x5fold', 0.05, LF_his, instdof, [-np.inf, np.inf, -np.inf, np.inf])
# # If do not need load-displacement diagram:
# # VisualFold(U_his[:, :interv:endicrm], truss, angles, 'none', 'miura5x5fold', 0.0001)

# Plot stored energy vs. pseudo time
# Red line is the total profile. Between red and cyan is the folding
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

