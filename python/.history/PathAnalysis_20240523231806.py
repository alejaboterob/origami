import numpy as np
from scipy.linalg import solve
from GlobalK_edu_ver import GlobalK_edu_ver
from nlsmgd import nlsmgd
from GetData import GetData


from scipy.io import loadmat

def PathAnalysis(truss, angles, F, b_lambda, MaxIcr):
    tol = 1e-6
    Maxitera = 50
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

    # F = F[0]
    while icrm < MaxIcr:
        icrm += 1
        itera = 0
        err = 1
        print(f'icrm = {icrm}, lambda = {lmd:6.4f}')
        while err > tol and itera < Maxitera:
            itera += 1
            IF, K = GlobalK_edu_ver(U, Node, truss, angles)

            # IF, K = GlobalK_edu_ver(U, Node, truss, angles)
            data = loadmat('data.mat')

            # The variables are stored in a dictionary
            # Access the variables using their MATLAB names
            IF = data['IF']
            K = data['K']


            R = lmd * F - IF.T
            MRS = np.column_stack((F, R.T))

     


            MUL[FreeDofs, :] = solve(K[FreeDofs, :][:, FreeDofs].A, MRS[FreeDofs, :])
            dUp = MUL[:, 0]
            dUr = MUL[:, 1]
            if itera == 1:
                dUr = np.zeros_like(dUr)
            dlmd = nlsmgd(icrm, itera, dUp, dUr, b_lambda)
            dUt = dlmd * dUp + dUr
            U = U + dUt
            err = np.linalg.norm(dUt[FreeDofs])
            lmd = lmd + dlmd
            print(f'    itera = {itera}, err = {err:6.4f}, dlambda = {dlmd:6.4f}')
            if err > 1e8:
                print('Divergence!')
                break

        if itera > 15:
            b_lambda = b_lambda / 2
            print('Reduce constraint radius!')
            icrm -= 1
            U = Uhis[:, max(icrm, 1) - 1]
            lmd = load_his[max(icrm, 1) - 1]
        elif itera < 3:
            print('Increase constraint radius!')
            b_lambda = b_lambda * 1.5
            Uhis[:, icrm - 1] = U
            load_his[icrm - 1] = lmd
            Exbari, FdAnglei, BdAnglei, LFdAnglei, LBdAnglei = GetData(U, Node, truss, angles)
            Data['Exbar'][:, icrm-1 ] = Exbari.T
            Data['FdAngle'][:, icrm - 1] = FdAnglei.T
            Data['LFdAngle'][:, icrm - 1] = LFdAnglei.T
            Data['BdAngle'][:, icrm - 1] = BdAnglei.T
            Data['LBdAngle'][:, icrm - 1] = LBdAnglei.T
        else:
            Uhis[:, icrm - 1] = U
            load_his[icrm - 1] = lmd
            Exbari, FdAnglei, BdAnglei, LFdAnglei, LBdAnglei = GetData(U, Node, truss, angles)
            Data['Exbar'][:, icrm - 1] = Exbari.T
            Data['FdAngle'][:, icrm - 1] = FdAnglei.T
            Data['LFdAngle'][:, icrm - 1] = LFdAnglei.T
            Data['BdAngle'][:, icrm - 1] = BdAnglei.T
            Data['LBdAngle'][:, icrm - 1] = LBdAnglei.T

    icrm += 1
    Uhis = Uhis[:, :icrm]
    load_his = load_his[:icrm]
    Data['Exbar'] = Data['Exbar'][:, :icrm]
    Data['FdAngle'] = Data['FdAngle'][:, :icrm]
    Data['LFdAngle'] = Data['LFdAngle'][:, :icrm]
    Data['BdAngle'] = Data['BdAngle'][:, :icrm]
    Data['LBdAngle'] = Data['LBdAngle'][:, :icrm]

    return Uhis, load_his, Data