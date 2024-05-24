import numpy as np

def icross(a, b):
    c = np.array([a[1,:]*b[2,:] - a[2,:]*b[1,:],
                 a[2,:]*b[0,:] - a[0,:]*b[2,:],
                 a[0,:]*b[1,:] - a[1,:]*b[0,:]])
    return c

