import numpy as np

def nlsmgd(step, ite, dup, dur, cmp, dupc1, numgsp, dupp1, sinal):
    """
    Modified generalized displacement control method.
    """
    if step == 1 and ite == 1:
        sinal = 1
        dupp1 = dup

    if ite == 1:
        sinal = np.sign(np.dot(dupp1, dup)) * sinal
        dupp1 = dup

    if ite == 1:
        if step == 1:
            dl = cmp
            numgsp = np.dot(dup, dup)
            dupc1 = dup
        else:
            gsp = numgsp / np.dot(dup, dup)
            dl = sinal * cmp * np.sqrt(gsp)
            dupc1 = dup
    else:
        dl = -np.dot(dupc1, dur) / np.dot(dupc1, dup)
    return dl, dupc1, numgsp, dupp1, sinal

