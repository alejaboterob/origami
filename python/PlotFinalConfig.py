import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, writers
# from PIL import Image
# from PlotOri import PlotOri
from plot_ori_panels import plot_panels, save_gif_PIL

def PlotFinalConfig(U_his, truss, angles, LF_his=None):
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


    mins = np.minimum(Node.min(axis=0), Node.min(axis=0))
    maxs = np.maximum(Node.max(axis=0), Node.max(axis=0))

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
    plot_panels(Nodew, Panel, ax=ax, **plot_kwargs)    
    plt.xticks([])
    plt.yticks([])
    ax.set_zticks([])
    ax.set(xticklabels=[], yticklabels=[], zticklabels=[])
    ax.auto_scale_xyz([mins[0], maxs[0]],
                    [mins[1], maxs[1]],
                    [mins[2], maxs[2]])    

    plt.show(block = False)

