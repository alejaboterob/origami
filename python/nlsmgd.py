import numpy as np

def nlsmgd(step, ite, dup, dur, cmp):
    '''The function `nlsmgd` implements the Modified Generalized Displacement Control Method for
    calculating displacement increments in a structural analysis.
    
    Parameters
    ----------
    step
        Step indicates the current step in the calculation process.
    ite
        The `ite` parameter in the `nlsmgd` function represents the iteration number. It is used to keep
    track of the current iteration within the function.
    dup
        The `dup` parameter in the `nlsmgd` function seems to represent a displacement vector. It is used
    to calculate various values and perform operations within the function.
    dur
        The parameter `dur` in the `nlsmgd` function represents the displacement vector. It is used in
    calculations to determine the displacement control value `dl` based on the input parameters and
    previous iterations within the function.
    cmp
        The `cmp` parameter in the `nlsmgd` function stands for the current displacement control parameter.
    It is used in the calculation of the displacement `dl` within the function.
    
    Returns
    -------
        The function `nlsmgd` is returning the variable `dl`.
    
    '''
    """
    Modified generalized displacement control method.
    """
    global dupp1, sinal
    
    if step == 1 and ite == 1:
        sinal = 1
        dupp1 = dup.copy()

    if ite == 1:
        sinal = np.sign(np.dot(dupp1, dup)) * sinal
        dupp1 = dup.copy()

    global dupc1, numgsp

    if ite == 1:
        if step == 1:
            dl = cmp
            numgsp = np.dot(dup, dup)
            dupc1 = dup.copy()
        else:
            gsp = numgsp / np.dot(dup, dup)
            dl = sinal * cmp * np.sqrt(gsp)
            dupc1 = dup.copy()
    else:
        dl = -np.dot(dupc1, dur) / np.dot(dupc1, dup)
    return dl

