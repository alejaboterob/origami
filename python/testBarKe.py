import numpy as np
from scipy.sparse import csr_matrix
from BarKe import BarKe
from Ogden import Ogden

import scipy.io


# Example test for the BarKe function
def test_BarKe():
    # Example inputs
    Ui_example = np.array([0.11449873, 0.60748948, 0.92399539, 0.02126101, 0.24440479,
       0.48081165])  # Adjust size and values based on expected input
    B_example = csr_matrix(np.array([[0.99377798, 0.51655026, 0.52905187, 0.89969584, 0.50178366,
        0.91173732]]))  # Adjust dimensions and density
    L_example = 10  # Adjust based on expected input
    E0 = 1e6

    CM_example = lambda Ex: Ogden(Ex, E0)  # Define bar material constitutive
    A_example = 1.0  # Cross-sectional area

    expected_Rbe = scipy.io.loadmat('Rbe.mat')['Rbe']
    expected_Kbe = scipy.io.loadmat('Kbe.mat')['Kbe']

    # Running the function
    _, Rbe, Kbe = BarKe(Ui_example, B_example, L_example, CM_example, A_example)
    print(Rbe-expected_Rbe.T)
    print(Kbe-expected_Kbe)

    # Check if the output matches the expected output
    assert np.allclose(Rbe, expected_Rbe.T), "Rbe does not match expected"
    assert np.allclose(Kbe, expected_Kbe), "Kbe does not match expected"

    print("BarKe function passed all tests")

# Run the test
test_BarKe()
