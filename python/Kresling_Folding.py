import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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

# Set boundary condition
indsupp = np.where(Node[:, 2] < 0.01)[0]
nsupp = len(indsupp)
Supp = np.zeros((nsupp + 2, 4))
Supp[0, 0] = indsupp[0]
Supp[0, 1:] = 1
Supp[1, 0] = indsupp[1]
Supp[1, 1:] = 1
Supp[2:, 0] = indsupp[2:]
Supp[2:, 1:3] = 1
Supp[2:, 3] = 1

m = Node.shape[0]
indp = np.where(np.abs(Node[:, 2] - np.max(Node[:, 2])) < 1e-5)[0]
npp = len(indp)
Load = np.zeros((npp, 4))
Load[:, 0] = indp
Load[:, 3] = -1

# Prepare data
truss, angles, F = PrepareData(Node, Panel, Supp, Load, BarMater, RotSpring, Kf, Kb, Abar)
truss.U0 = np.zeros(3 * Node.shape[0])
U_his, LF_his, Data = PathAnalysis(truss, angles, F, blam, MaxIcr)
U_his = np.real(U_his)
LF_his = np.real(LF_his)

# Visualize simulation
instdof = -indp[0] * 3
interv = 1
endicrm = U_his.shape[1]
VisualFold(U_his[:, ::interv], truss, angles, 'none', 'kreslingfold', 0.0001, LF_his, instdof, [-np.inf, np.inf, -np.inf, np.inf])

# Plot stored energy vs. pseudo time
STAT = PostProcess(Data, truss, angles)
plt.figure()
plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.PE.strain, 'r-', linewidth=2)
plt.grid()
plt.hold(True)
plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.bend.PE + STAT.bar.PE, 'c-')
plt.plot(np.arange(1, U_his.shape[1] + 1), STAT.bar.PE, 'm-')
plt.xlabel('Increment Number (Pseudo-time)', fontsize=14)
plt.ylabel('Stored Energy', fontsize=14)

# Plot final configuration
Ux = U_his[:, -1]
Nodew = Node.copy()
Nodew[:, 0] += Ux[::3]
Nodew[:, 1] += Ux[1::3]
Nodew[:, 2] += Ux[2::3]
plt.figure()
PlotOri(Nodew, angles.Panel, truss.Trigl)
plt.axis('equal')
plt.axis('off')
plt.camproj('perspective')
plt.light()
plt.view(117, 18)
plt.rotate3d()

