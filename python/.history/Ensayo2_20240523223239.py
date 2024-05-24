import numpy as np
from ConfigMiura import ConfigMiura
from Ogden import Ogden
from PrepareData import PrepareData
from PathAnalysis import PathAnalysis
from EnhancedLinear import EnhancedLinear

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



from scipy.io import loadmat
import numpy as np
from scipy.linalg import solve
from GlobalK_edu_ver import GlobalK_edu_ver
from nlsmgd import nlsmgd
from GetData import GetData

data = loadmat('data.mat')

# The variables are stored in a dictionary
# Access the variables using their MATLAB names
IF2 = data['IF']
K2 = data['K']

print("IF:", IF2)
print("K:", K2)

b_lambda = blam

tol = 1e-6
Maxite = 50
Node = truss['Node']
AllDofs = np.arange(1, 3 * Node.shape[0])  # sin +1 igual a matlab 
U = truss['U0']
Uhis = np.zeros((3 * Node.shape[0], MaxIcr))
Data = {}
Data['Exbar'] = np.zeros((truss['Bars'].shape[0], MaxIcr))
Data['FdAngle'] = np.zeros((angles['fold'].shape[0], MaxIcr))
Data['LFdAngle'] = Data['FdAngle'].copy()
Data['BdAngle'] = np.zeros((angles['bend'].shape[0], MaxIcr))
Data['LBdAngle'] = Data['BdAngle'].copy()

FreeDofs = np.setdiff1d(AllDofs, truss['FixedDofs'])
lmd = 0
icrm = 0
MUL = np.column_stack((U, U))
load_his = np.zeros(MaxIcr)


R = lmd * F - IF2.T
MRS = np.column_stack((F, R.T))

ite==1

MUL[FreeDofs, :] = solve(K2[FreeDofs, :][:, FreeDofs].A, MRS[FreeDofs, :])
dUp = MUL[:, 0]
dUr = MUL[:, 1]
if ite == 1:
    dUr = np.zeros_like(dUr)
dlmd = nlsmgd(icrm, ite, dUp, dUr, blam)
dUt = dlmd * dUp + dUr
U = U + dUt
err = np.linalg.norm(dUt[FreeDofs])
lmd = lmd + dlmd
print(f'    ite = {ite}, err = {err:6.4f}, dlambda = {dlmd:6.4f}')


if ite > 15:
    b_lambda = b_lambda / 2
    print('Reduce constraint radius!')
    icrm -= 1
    U = Uhis[:, max(icrm, 1) - 1]
    lmd = load_his[max(icrm, 1) - 1]
elif ite < 3:
    print('Increase constraint radius!')
    b_lambda = b_lambda * 1.5
    Uhis[:, icrm - 1] = U
    load_his[icrm - 1] = lmd
    Exbari, FdAnglei, BdAnglei, LFdAnglei, LBdAnglei = GetData(U, Node, truss, angles)
    Data['Exbar'][:, icrm - 1] = Exbari
    Data['FdAngle'][:, icrm - 1] = FdAnglei
    Data['LFdAngle'][:, icrm - 1] = LFdAnglei
    Data['BdAngle'][:, icrm - 1] = BdAnglei
    Data['LBdAngle'][:, icrm - 1] = LBdAnglei
else:
    Uhis[:, icrm - 1] = U
    load_his[icrm - 1] = lmd
    Exbari, FdAnglei, BdAnglei, LFdAnglei, LBdAnglei = GetData(U, Node, truss, angles)
    Data['Exbar'][:, icrm - 1] = Exbari
    Data['FdAngle'][:, icrm - 1] = FdAnglei
    Data['LFdAngle'][:, icrm - 1] = LFdAnglei
    Data['BdAngle'][:, icrm - 1] = BdAnglei
    Data['LBdAngle'][:, icrm - 1] = LBdAnglei

icrm += 1
Uhis = Uhis[:, :icrm]
load_his = load_his[:icrm]
Data['Exbar'] = Data['Exbar'][:, :icrm]
Data['FdAngle'] = Data['FdAngle'][:, :icrm]
Data['LFdAngle'] = Data['LFdAngle'][:, :icrm]
Data['BdAngle'] = Data['BdAngle'][:, :icrm]
Data['LBdAngle'] = Data['LBdAngle'][:, :icrm]