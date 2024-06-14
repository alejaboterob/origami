import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, writers
# from PIL import Image
# from PlotOri import PlotOri
from plot_ori_panels import plot_panels, save_gif_PIL

def VisualFold(U_his, truss, angles, LF_his=None):
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

    plot_kwargs_wire = {"linewidths": 1, "edgecolors": "#3c3c3c22",
               "facecolors": "#d86a96", "alpha": 0.0, "zorder": 4,
               "zsort": "max"}
    plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
               "facecolors": "#d86a96", "zorder": 5,
               "zsort": "min"}
    



    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_kwargs_wire = {"linewidths": 0.5, "edgecolors": "#858585",
                   "facecolors": "#d86a96", "alpha": 0.2, "zorder": 0}
    plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
                   "facecolors": "#d86a96", "zorder": 10.5}
    
    files = []
    Node_history = []

    for i in range(U_his.shape[1]):
        plt.cla()
        U = U_his[:,i]
        # plt.view(35, 30)
        Nodew = Node.copy()
        Nodew[:, 0] = Node[:, 0] + U[::3]
        Nodew[:, 1] = Node[:, 1] + U[1::3]
        Nodew[:, 2] = Node[:, 2] + U[2::3]
        plot_panels(Node, Panel, ax=ax, **plot_kwargs_wire)
        plot_panels(Nodew, Panel, ax=ax, **plot_kwargs)    
        # plt.xticks([0, 5, 10, 15, 20])
        # plt.yticks([0, 5, 10, 15, 20])
        # ax.set_zticks([0, 2])
        plt.xticks([])
        plt.yticks([])
        ax.set_zticks([])
        ax.set(xticklabels=[], yticklabels=[], zticklabels=[])
        ax.auto_scale_xyz([mins[0], maxs[0]],
                        [mins[1], maxs[1]],
                        [mins[2], maxs[2]])
        plt.axis("image")
        file = f"ori_{str(i).zfill(2)}.png"
        plt.savefig(file)
        files.append(file)
        Node_history.append(Nodew)
    
    data = {'Node': Node, 'Node_history': Node_history, 'Panel': Panel}
    np.save('data_completa.npy', data)

    plt.show()

    save_gif_PIL("ori_anim.gif", files, fps=5, loop=0)
    # [os.remove(file) for file in files]
