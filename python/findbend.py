import numpy as np

def findbend(Panel, Node):
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