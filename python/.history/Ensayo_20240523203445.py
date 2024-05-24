# MIURA FOLDING
import numpy as np

# Define geometry and material parameters
x_divs = 1  # Number of unit cells in horizontal direction
y_divs = 1  # Number of unit cells in vertical direction
# Geometry of the Miura: a, b are the edge lengths of each parallelogram
# panel; fdang controls the folding angle; theta is the panel angle
theta = 60
a = 2
b = 2
gmma = 15

# Maximum increment number & initial load factor
MaxIcr = 60
blam = 0.5

# Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
# E0-stretching stiffness; Abar-bar area (assume uniform)
Kf = 1e-1
Kb = Kf * 1e5
E0 = 1e6
Abar = 1e-1

# Left and right limits for the linear range of rotational stiffness
limlft = 0.1
limrht = 360 - 0.1    


theta = np.pi * theta / 180
gmma = gmma * np.pi / 180
numx = 2 * x_divs
numy = 2 * y_divs
h = a * np.sin(theta) * np.sin(gmma)
s = b * np.cos(gmma) * np.tan(theta) / np.sqrt(1 + np.cos(gmma) ** 2 * np.tan(theta) ** 2)
l = a * np.sqrt(1 - np.sin(gmma) ** 2 * np.sin(theta) ** 2)
v = b / np.sqrt(1 + np.cos(gmma) ** 2 * np.tan(theta) ** 2)
X = np.repeat(s * np.arange(numx + 1), numy + 1).reshape(numy + 1, numx + 1)
Y = np.repeat(l * np.arange(numy + 1), numx + 1).reshape(numy + 1, numx + 1)
Y[:, 1::2] = Y[:, 1::2] + v
Z = np.zeros_like(X)
Z[1::2, :] = h

# Reshape each array to a column vector
X_reshaped = X.T.reshape(-1, 1)
Y_reshaped = Y.T.reshape(-1, 1)
Z_reshaped = Z.T.reshape(-1, 1)

# Concatenate the reshaped arrays column-wise
NODE = np.hstack((X_reshaped, Y_reshaped, Z_reshaped))

print(Y.flatten())	
print(NODE.shape)
print(Z)

from scipy.io import savemat

# Assuming you have already calculated the variables NODE, PANEL, and BDRY
# and they are numpy arrays
variables = {
    'X2': X_reshaped,
    'Y2': Y_reshaped,
    'Z2': Z_reshaped,
    'NODE2': NODE,
    # 'PANEL2': Panel,
    # 'BDRY2': BDRY
}

# Save the variables to a .mat file
savemat('variables.mat', variables)



# variables = {
#     'X2': X,
#     'Y2': Y,
# }
# savemat('variables.mat', variables)


# np.savetxt('B.csv', B.todense(), delimiter = ',') 
