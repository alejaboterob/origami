import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from GetDiSym import GetDiSym  
from Ogden import Ogden
from EnhancedLinear import EnhancedLinear
from PrepareData import PrepareData
from PathAnalysis import PathAnalysis   
from VisualFold import VisualFold   
from PostProcess import PostProcess 
from PlotOri import PlotOri 
from PlotFinalConfig import PlotFinalConfig
from displacement import displacement
from GraphPostProcess import GraphPostProcess

# Define geometry and material
N = 8
h = 1
lyr = 6
phi = 2 * np.pi / 8

MaxIcr = 180
blam = 0.032

Kf = 1e-3
Kb = Kf
E0 = 5e3
Abar = 0.1

limlft = 45
limrht = 315

# Left and right limits for the linear range of rotational stiffness
[Node, Panel] = GetDiSym(N,h,lyr,phi)


# Define bar material constitutive and rotational spring functions
BarMater = lambda Ex: Ogden(Ex, E0)  # Define bar material constitutive
RotSpring = lambda he, h0, kpi, L0: EnhancedLinear(he, h0, kpi, L0, limlft, limrht)

# Set boundary condition
indsupp = np.where(Node[:, 2] < 0.01)[0]
nsupp = len(indsupp)
Supp = np.ones((nsupp, 4))
Supp[0:, 0] = indsupp[0:]

m = Node.shape[0]
indp = np.where(np.abs(Node[:, 2] - np.max(Node[:, 2])) < 1e-5)[0]
npp = len(indp)
Load = np.zeros((npp, 4))
Load[:, 0] = indp
Load[:, 3] = -1

# Prepare data

truss, angles, F = PrepareData(Node, Panel, Supp, Load, BarMater, RotSpring, Kf, Kb, Abar)
truss['U0'] = np.zeros(3 * Node.shape[0])
U_his, LF_his, Data = PathAnalysis(truss, angles, F, blam, MaxIcr)
U_his = np.real(U_his)
LF_his = np.real(LF_his)

# Visualize simulation
instdof = -indp[0] * 3
interv = 1
endicrm = U_his.shape[1]
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


# # Plot stored energy vs. pseudo time
# STAT = PostProcess(Data, truss, angles)
# plt.figure()
# plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.PE.strain, 'r-', linewidth=2)
# plt.grid()
# plt.hold(True)
# plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.bend.PE + STAT.bar.PE, 'c-')
# plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.bar.PE, 'm-')
# plt.xlabel('Increment Number (Pseudo-time)', fontsize=14)
# plt.ylabel('Stored Energy', fontsize=14)

# # Plot final configuration
# Ux = U_his[:, -1]
# Nodew = Node.copy()
# Nodew[:, 0] += Ux[::3]
# Nodew[:, 1] += Ux[1::3]
# Nodew[:, 2] += Ux[2::3]
# plt.figure()
# PlotOri(Nodew, angles.Panel, truss.Trigl)
# plt.axis('equal')
# plt.axis('off')
# plt.camproj('perspective')
# plt.light()
# plt.view(117, 18)
# plt.rotate3d()

