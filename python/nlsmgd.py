import numpy as np

def nlsmgd(step, ite, dup, dur, cmp):
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

