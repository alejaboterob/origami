import numpy as np

def PostProcess(Data, truss, angles):
    '''The function `PostProcess` calculates various statistics related to truss and angles based on input
    data.
    
    Parameters
    ----------
    Data
        The `PostProcess` function takes in three parameters: `Data`, `truss`, and `angles`. Here is a
    brief description of each parameter:
    truss
        The `truss` parameter seems to be a dictionary containing information related to the truss
    structure. It likely includes keys such as 'CM', 'L', 'A', and possibly others that are used within
    the `PostProcess` function for calculations.
    angles
        The `angles` parameter in the `PostProcess` function likely contains information related to the
    angles of the truss structure being analyzed. It seems to include keys such as 'pf0', 'kpf', 'pb0',
    'kpb', 'bend', and 'fold'. These keys
    
    Returns
    -------
        The function `PostProcess` returns a dictionary `STAT` containing the following keys and values:
    
    '''
    Sx_bar, _, Wb = truss['CM'](Data['Exbar'])
    
    Rspr_fd = np.zeros(Data['FdAngle'].shape)
    Efold = Rspr_fd.copy()
    
    Rspr_bd = np.zeros(Data['BdAngle'].shape)
    Ebend = Rspr_bd.copy()
    
    for i in range(Data['FdAngle'].shape[1]):
        Rspr_fdi, _, Efoldi = angles['CM'](Data['FdAngle'][:, i], angles['pf0'], angles['kpf'], truss['L'][angles['bend'].shape[0]:angles['bend'].shape[0] + angles['fold'].shape[0]])
        Rspr_bdi, _, Ebendi = angles['CM'](Data['BdAngle'][:, i], angles['pb0'], angles['kpb'], truss['L'][:angles['bend'].shape[0]])
        Rspr_fd[:, i] = Rspr_fdi
        Efold[:, i] = Efoldi
        Rspr_bd[:, i] = Rspr_bdi
        Ebend[:, i] = Ebendi
    
    STAT = {}
    STAT['bar'] = {}
    STAT['bar']['Sx'] = Sx_bar
    STAT['bar']['W'] = np.diag(truss['L'] * truss['A']) @ Wb
    STAT['bar']['PE'] = np.sum(STAT['bar']['W'], axis=0)
    
    STAT['fold'] = {}
    STAT['fold']['RM'] = Rspr_fd
    STAT['fold']['E'] = Efold
    STAT['fold']['PE'] = np.sum(STAT['fold']['E'], axis=0)
    
    STAT['bend'] = {}
    STAT['bend']['RM'] = Rspr_bd
    STAT['bend']['E'] = Ebend
    STAT['bend']['PE'] = np.sum(STAT['bend']['E'], axis=0)
    
    STAT['PE'] = {}
    STAT['PE']['strain'] = STAT['bar']['PE'] + STAT['fold']['PE'] + STAT['bend']['PE']
    
    return STAT

