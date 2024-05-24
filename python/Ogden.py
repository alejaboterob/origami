import numpy as np
from scipy.sparse import issparse

def Ogden(Ex, C0):
    # Ogden hyperelastic constitutive model for bar elements
    alfa = [5, 1]  # Specify parameters  
    pstr = np.real(np.sqrt(2*Ex+1))
    C0 = (pstr<1)*1*C0+(~(pstr<1))*1*C0
    Et = C0/(alfa[0]-alfa[1])*((alfa[0]-2)*pstr**(alfa[0]-4)-(alfa[1]-2)*pstr**(alfa[1]-4))
    Sx = C0/(alfa[0]-alfa[1])*(pstr**(alfa[0]-2)-pstr**(alfa[1]-2))
    Wb = C0/(alfa[0]-alfa[1])*((pstr**alfa[0]-1)/alfa[0]-(pstr**alfa[1]-1)/alfa[1])
    return Sx, Et, Wb


