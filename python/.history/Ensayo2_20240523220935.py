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


R = lmd * F - IF.T
MRS = np.column_stack((F, R.T))


MUL[FreeDofs, :] = solve(K2[FreeDofs, :][:, FreeDofs].A, MRS[FreeDofs, :])
dUp = MUL[:, 0]
dUr = MUL[:, 1]
if iter == 1:
    dUr = np.zeros_like(dUr)
dlmd = nlsmgd(icrm, iter, dUp, dUr, b_lambda)
dUt = dlmd * dUp + dUr
U = U + dUt
err = np.linalg.norm(dUt[FreeDofs])
lmd = lmd + dlmd
print(f'    iter = {iter}, err = {err:6.4f}, dlambda = {dlmd:6.4f}')


if iter > 15:
    b_lambda = b_lambda / 2
    print('Reduce constraint radius!')
    icrm -= 1
    U = Uhis[:, max(icrm, 1) - 1]
    lmd = load_his[max(icrm, 1) - 1]
elif iter < 3:
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