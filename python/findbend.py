import numpy as np

def findbend(Panel, Node):
    '''The `findbend` function in Python calculates the bending points of panels based on the distances
    between nodes.
    
    Parameters
    ----------
    Panel
        The `Panel` parameter seems to be a list of lists where each inner list represents a panel with
    four nodes. Each inner list contains the indices of the four nodes that make up the panel.
    Node
        The `Node` parameter in the `findbend` function likely represents an array containing the
    coordinates of nodes in a structural panel or mesh. These coordinates are used to calculate the
    lengths of edges in the panel.
    
    Returns
    -------
        The function `findbend` returns a numpy array containing the indices of nodes that define the
    bending sequence for each panel in the input `Panel` array. The returned array has a shape of (n, 4)
    where n is the number of panels with valid bending sequences.
    
    '''
    bend = np.zeros((len(Panel), 4), dtype=int)
    for i in range(len(Panel)):
        if len(Panel[i]) == 4:
            L1 = np.linalg.norm(Node[Panel[i][0]] - Node[Panel[i][2]])   #cada uno con -1 igual a matlab
            L2 = np.linalg.norm(Node[Panel[i][3]] - Node[Panel[i][1]])   #cada uno con -1 igual a matlab
            if L1 > L2:
                lclbend = [Panel[i][1], Panel[i][3], Panel[i][0], Panel[i][2]]
            else:
                lclbend = [Panel[i][0], Panel[i][2], Panel[i][1], Panel[i][3]]
            bend[i, :] = lclbend
    bend = bend[np.sum(bend, axis=1) != 0, :]
    return bend