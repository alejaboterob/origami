import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def PlotOri(node, panel, trigl, lstyle='-', alfa=1, color='g'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Default settings
    lstl = lstyle
    al = alfa
    cc = color
    edge_color = (1 - al) * np.array([1, 1, 1])
    
    # Function to create patches from indices
    def create_patches(indices):
        patches = []
        for index_set in indices:
            polygon = node[index_set - 1]  # -1 for zero index
            patches.append(polygon)
        return patches

    # Plot triangles if trigl is not empty
    if trigl.size > 0:
        trig_patches = create_patches(trigl)
        pc = Poly3DCollection(trig_patches, facecolor=cc, linestyle='none', 
                             edgecolor=edge_color, linewidth=1)
        ax.add_collection3d(pc)

    # Organize panels by size
    panelsize = np.array([len(p) for p in panel])
    ptri = [panel[i] for i in np.where(panelsize == 3)[0]]
    pquad = [panel[i] for i in np.where(panelsize == 4)[0]]

    # Plot panels based on the trigl array
    if trigl.size == 0:
        ptri_patches = create_patches(ptri)
        pquad_patches = create_patches(pquad)
        pc_tri = Poly3DCollection(ptri_patches, facecolor=cc, linestyle=lstl,
                                 linewidth=1, edgecolor=edge_color)
        pc_quad = Poly3DCollection(pquad_patches, facecolor=cc, linestyle=lstl,
                                  linewidth=1, edgecolor=edge_color)
        ax.add_collection3d(pc_tri)
        ax.add_collection3d(pc_quad)
    else:
        ptri_patches = create_patches(ptri)
        pquad_patches = create_patches(pquad)
        pc_tri = Poly3DCollection(ptri_patches, facecolor='none', linestyle=lstl,
                                 linewidth=1, edgecolor=edge_color)
        pc_quad = Poly3DCollection(pquad_patches, facecolor='none', linestyle=lstl,
                                  linewidth=1, edgecolor=edge_color)
        ax.add_collection3d(pc_tri)
        ax.add_collection3d(pc_quad)

    ax.set_xlim([node[:, 0].min(), node[:, 0].max()])
    ax.set_ylim([node[:, 1].min(), node[:, 1].max()])
    ax.set_zlim([node[:, 2].min(), node[:, 2].max()])
    plt.show()
