import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, writers
# from PIL import Image
# from PlotOri import PlotOri
from plot_ori_panels import plot_panels, save_gif_PIL

def PlotFinalConfig(U_his, truss, angles, LF_his=None):
    '''The function `PlotFinalConfig` plots the final configuration of a truss structure with optional
    features for recording and displaying load vs. displacement diagrams.
    
    Parameters
    ----------
    U_his
        It seems like the code snippet you provided is a function named `PlotFinalConfig` that is used to
    plot a truss structure with some specified configurations. The function takes several parameters
    including `U_his`, `truss`, `angles`, and optionally `LF_his`.
    truss
        The `truss` parameter in the `PlotFinalConfig` function seems to represent the truss structure in
    the simulation. It likely contains information about the nodes and elements of the truss. The `Node`
    and `Trigl` keys within the `truss` dictionary probably store the
    angles
        The `angles` parameter in the `PlotFinalConfig` function seems to represent the angles of the truss
    panels. It is used to plot the panels in the 3D visualization. The `Panel` variable extracted from
    the `angles` dictionary likely contains information about the panels such as their vertices
    LF_his
        LF_his is a NumPy array that represents the load history of the truss structure. If provided, it is
    used to calculate the total load at each time step for plotting load vs. displacement diagram.
    
    '''
    """
    Record the simulation to file if needed:
    recordtype = 'none': do not save the simulation
    recordtype = 'video': save simulation in MP4 format;
    recordtype = 'imggif': save simulatrion in GIF format.
    pausement: pause time between each frame (in seconds).
               If recordtype = 'none': use a small number such as 0.0001;
               Otherwise, pausetime = 1/fps;
    If input data does not include 'load_his', 'instdof', 'axislim', the
    function does not plot load vs. displacement diagram.
    axislim: Axis limits (bounding box) for load vs. displacement diagram.
             format: [xmin,xmax,ymin,ymax].
    instdof: specify the DOF of interest for displacement measure.
    """

    Node = truss['Node']
    Trigl = truss['Trigl']
    Panel = angles['Panel']

    if LF_his is not None and len(LF_his.shape) > 1:
        LF_his = np.sum(LF_his, axis=1)




    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_kwargs_wire = {"linewidths": 0.5, "edgecolors": "#858585",
                   "facecolors": "#d86a96", "alpha": 0.2, "zorder": 4}
    plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
                   "facecolors": "#d86a96", "zorder": 5}

    Ux = U_his[:, -1]

    Nodew = Node.copy()
    Nodew[:, 0] = Node[:, 0] + Ux[0::3]
    Nodew[:, 1] = Node[:, 1] + Ux[1::3]
    Nodew[:, 2] = Node[:, 2] + Ux[2::3]

    mins = np.minimum(Nodew.min(axis=0), Nodew.min(axis=0))
    maxs = np.maximum(Nodew.max(axis=0), Nodew.max(axis=0))

    plot_panels(Nodew, Panel, ax=ax, **plot_kwargs)    
    plt.xticks([])
    plt.yticks([])
    ax.auto_scale_xyz([mins[0], maxs[0]],
                [mins[1], maxs[1]],
                [mins[2], maxs[2]])    
    ax.set_zticks([])
    ax.axis('equal')
    ax.set(xticklabels=[], yticklabels=[], zticklabels=[])


    plt.show(block = False)

