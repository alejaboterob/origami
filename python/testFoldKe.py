import numpy as np
from scipy.sparse import csr_matrix
from FoldKe import FoldKe  # Ensure you have this import if FoldKe is a separate module
import scipy.io

def test_FoldKe():
    # Example inputs
    Nodenw_example = np.array([...])  # Placeholder
    bend_example = np.array([...])  # Placeholder
    angles_example = {'kpb': np.array([...]), 'pb0': np.array([...]), 'CM': ...}  # Placeholder, adjust accordingly
    Lbend_example = np.array([...])  # Placeholder
    d_el = 0  # Assuming 'd_el' is an index, adjust based on actual usage

    # Load expected outputs from MATLAB .mat files
    expected_Rpe = scipy.io.loadmat('Rpe.mat')['Rpe']
    expected_Kpe = scipy.io.loadmat('Kpe.mat')['Kpe']

    # Running the function
    _, Rpe, Kpe = FoldKe(Nodenw_example, bend_example, angles_example['kpb'], angles_example['pb0'][d_el], Lbend_example[d_el], angles_example['CM'])
    
    print("Difference in Rpe:", Rpe - expected_Rpe.T)
    print("Difference in Kpe:", Kpe - expected_Kpe)

    # Check if the output matches the expected output
    assert np.allclose(Rpe, expected_Rpe.T), "Rpe does not match expected"
    assert np.allclose(Kpe, expected_Kpe), "Kpe does not match expected"

    print("FoldKe function passed all tests")

# Run the test
test_FoldKe()
