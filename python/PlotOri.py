import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def PlotOri(node, panel, trigl, lstyle='-', alfa=1, color='g'):
    '''The function `PlotOri` in Python is used to plot 3D triangles and panels based on input node, panel,
    and trigl data with customizable line style, transparency, and color.
    
    Parameters
    ----------
    node
        The `node` parameter in the `PlotOri` function seems to represent the coordinates of nodes in a 3D
    space. These nodes likely define the vertices of the geometric shapes (triangles and panels) that
    are being plotted in the function. The `node` parameter is used to create
    panel
        The `panel` parameter in the `PlotOri` function seems to represent a list of panels in a 3D space.
    Each panel is defined by a set of points that form a polygon. The panels are organized based on
    their size, with triangles and quadrilaterals separated into different lists
    trigl
        The `trigl` parameter in the `PlotOri` function seems to represent an array containing indices of
    triangles. These triangles are used to create 3D patches that will be plotted in the visualization.
    The function checks if `trigl` is not empty and then creates patches for
    lstyle, optional
        The `lstyle` parameter in the `PlotOri` function is used to specify the line style for plotting the
    panels. The default value is `'-'`, which represents a solid line style. You can change this
    parameter to customize the appearance of the lines in the plot. Some common line styles
    alfa, optional
        The `alfa` parameter in the `PlotOri` function seems to represent the transparency or alpha value
    of the plotted elements. A value of 1 indicates full opacity, while a value closer to 0 indicates
    increasing transparency.
    color, optional
        The `color` parameter in the `PlotOri` function is used to specify the color of the plotted
    polygons. The default color is set to 'g', which represents green. You can change this parameter to
    any valid color string in Matplotlib to customize the color of the plotted polygons. Some
    
    '''
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
            polygon = node[index_set]  # -1 for zero index
            patches.append(polygon)
        return patches

    # Plot triangles if trigl is not empty
    if trigl.size > 0:
        trig_patches = create_patches(trigl)
        pc = Poly3DCollection(trig_patches, facecolor=cc, linestyle='-', 
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
