import numpy as np
from scipy.linalg import solve
from GlobalK_edu_ver import GlobalK_edu_ver
from nlsmgd import nlsmgd
from GetData import GetData

from GlobalK_fast_ver import GlobalK_fast_ver
from scipy.io import loadmat

def PathAnalysis(truss, angles, F, b_lambda, MaxIcr):
    '''The function `PathAnalysis` performs a path-following analysis on a truss structure under specified
    loading conditions and constraints, storing the results at each iteration.
    
    Parameters
    ----------
    truss
        The `truss` parameter in the `PathAnalysis` function seems to be a dictionary containing
    information about the truss structure. It likely includes keys such as 'Node', 'Bars', 'FixedDofs',
    and 'U0', which are used within the function for calculations related to the tr
    angles
        The `angles` parameter seems to contain information related to the angles of the truss structure.
    It likely includes data about fold angles and bend angles. The function `PathAnalysis` uses this
    information along with other parameters to perform some calculations related to the truss structure.
    If you have specific questions or
    F
        The `F` parameter in the `PathAnalysis` function represents the applied load on the truss
    structure. It is a vector that contains the magnitudes of the applied loads on each degree of
    freedom of the truss.
    b_lambda
        The `b_lambda` parameter in the `PathAnalysis` function seems to be related to the constraint
    radius used in the analysis. It is adjusted during the iterations based on certain conditions within
    the function. The value of `b_lambda` affects the convergence behavior of the algorithm and is
    modified dynamically within the
    MaxIcr
        The `MaxIcr` parameter in the `PathAnalysis` function represents the maximum number of iterations
    allowed for the path analysis algorithm. This parameter controls how many iterations the algorithm
    will perform before terminating and returning the results.
    
    Returns
    -------
        The function `PathAnalysis` returns three main outputs:
    1. `Uhis`: A matrix containing the displacement history of the truss structure at each iteration.
    2. `load_his`: An array containing the load history of the truss structure at each iteration.
    3. `Data`: A dictionary containing various data arrays related to the truss structure, such as
    'Exbar', 'Fd
    
    '''
    tol = 1e-6
    Maxitera = 50
    Node = truss['Node']
    AllDofs = np.arange(0, 3 * Node.shape[0])  # sin +1 igual a matlab 
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
            # IF, K = GlobalK_fast_ver(U, Node, truss, angles)


            R = lmd * F - IF.T
            R[np.isnan(R)] = 0
            MRS = np.column_stack((F, R.T))

            MUL[FreeDofs, :] = solve(K[np.ix_(FreeDofs, FreeDofs)].A, MRS[FreeDofs, :])
            dUp = MUL[:, 0]
            dUr = MUL[:, 1]
            if itera == 1:
                dUr = np.zeros_like(dUr)
                # dupc12 = []   
                # numgsp = []
                # dupp1 = []
                # sinal = []

            # dlmd, dupc12, numgsp, dupp1, sinal = nlsmgd(icrm, itera, dUp, dUr, b_lambda,dupc12, numgsp, dupp1, sinal)
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