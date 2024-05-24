import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

def PlotOri(node, panel, trigl, lstyle='-', alfa=1, color='g'):
    fig, ax = plt.subplots()
    
    # Default settings
    lstl = lstyle
    al = alfa
    cc = color
    edge_color = (1 - al) * np.array([1, 1, 1])

    # Function to create patches from indices
    def create_patches(indices):
        patches = []
        for index_set in indices:
            polygon = Polygon(node[index_set - 1], closed=True)  # -1 for zero index
            patches.append(polygon)
        return patches

    # Plot triangles if trigl is not empty
    if trigl.size > 0:
        trig_patches = create_patches(trigl)
        pc = PatchCollection(trig_patches, facecolor=cc, linestyle='none', 
                             edgecolor=edge_color, linewidth=1)
        ax.add_collection(pc)

    # Organize panels by size
    panelsize = np.array([len(p) for p in panel])
    ptri = [panel[i] for i in np.where(panelsize == 3)[0]]
    pquad = [panel[i] for i in np.where(panelsize == 4)[0]]

    # Plot panels if trigl is empty
    if trigl.size == 0:
        ptri_patches = create_patches(ptri)
        pquad_patches = create_patches(pquad)
        pc_tri = PatchCollection(ptri_patches, facecolor=cc, linestyle=lstl,
                                 linewidth=1, edgecolor=edge_color)
        pc_quad = PatchCollection(pquad_patches, facecolor=cc, linestyle=lstl,
                                  linewidth=1, edgecolor=edge_color)
        ax.add_collection(pc_tri)
        ax.add_collection(pc_quad)
    else:
        ptri_patches = create_patches(ptri)
        pquad_patches = create_patches(pquad)
        pc_tri = PatchCollection(ptri_patches, facecolor='none', linestyle=lstl,
                                 linewidth=1, edgecolor=edge_color)
        pc_quad = PatchCollection(pquad_patches, facecolor='none', linestyle=lstl,
                                  linewidth=1, edgecolor=edge_color)
        ax.add_collection(pc_tri)
        ax.add_collection(pc_quad)

    plt.axis('equal')
    plt.show()



# import numpy as np
# import matplotlib.pyplot as plt

# def PlotOri(Node, Panel, Trigl, lstyle='-', alfa=1, color='g'):
#     if len(Node) < 4:
#         lstl = '-'
#         al = 1
#         cc = 'g'
#     else:
#         lstl = lstyle
#         al = alfa
#         cc = color
    
#     # if Trigl.any():
#         # plt.patch(faces=Trigl, vertices=Node, facecolor=cc, linestyle='none', facelighting='flat', edgecolor=(1-al)*np.ones(3))
    
#     # plt.hold(True)
    
#     Panelsize = [len(p) for p in Panel]
#     Ptri = [Panel[i] for i in np.where(np.array(Panelsize)==3)[0]]
#     Pquad = [Panel[i] for i in np.where(np.array(Panelsize)==4)[0]]
    
#     # if Trigl.any():
#         # plt.patch(faces=np.concatenate(Ptri), vertices=Node, facecolor='none', linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
#         # plt.patch(faces=np.concatenate(Pquad), vertices=Node, facecolor='none', linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
#     # else:
#         # plt.patch(faces=np.concatenate(Ptri), vertices=Node, facecolor=cc, linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
#         # plt.patch(faces=np.concatenate(Pquad), vertices=Node, facecolor=cc, linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
    
#     # plt.hold(False)

