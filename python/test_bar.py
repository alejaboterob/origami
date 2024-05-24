import numpy as np
import matplotlib.pyplot as plt
from Ogden import Ogden

E0 = 1
pstr = np.arange(0.1, 1.91, 0.01)

# Neo-Hookean
# alfa = [2, 0]

# St. Venant
# alfa = [4, 2] 

# Others
# alfa = [5, 1]

Ex = 0.5*(pstr**2 - 1)
Sx, Et, Wb = Ogden(Ex, E0)

fig = plt.figure(21)
plt.plot(pstr, Sx, linewidth=2)
plt.xlabel(r'$\lambda$', fontsize=14)
plt.ylabel(r'$S_{xx}$', fontsize=14)
# plt.axis([0.7, 1.3, -np.inf, np.inf])
plt.show()

from scipy.io import savemat
variables = {
    'pstr2': pstr,
    'Sx2': Sx,
    'Et2': Et,
    'Wb2': Wb,
}
# Save the variables to a .mat file
savemat('variables.mat', variables)