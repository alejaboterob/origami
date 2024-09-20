# -*- coding: utf-8 -*-
"""
Visualize origami panels in Matplotlib
@author: Nicolás Guarín-Zapata
@date: June 2024
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import LightSource
from PIL import Image
import os


def save_gif_PIL(outfile, files, fps=5, loop=0):
    '''The function `save_gif_PIL` saves a list of PNG files as a GIF with specified frames per second and
    looping options using the Python Imaging Library (PIL).
    
    Parameters
    ----------
    outfile
        The `outfile` parameter is a string that represents the path to the output file where the GIF will
    be saved.
    files
        The `files` parameter in the `save_gif_PIL` function is a list of paths with the PNG files that you
    want to include in the GIF. You should provide a list of file paths as strings to this parameter
    when calling the function.
    fps, optional
        The `fps` parameter in the `save_gif_PIL` function stands for frames per second. It determines how
    many frames are displayed per second in the resulting GIF animation.
    loop, optional
        The `loop` parameter in the `save_gif_PIL` function specifies the number of times the GIF should
    loop.
    
    '''
    """Helper function for saving GIFs
    
    Parameters
    ----------
    outfile : string
        Path to the output file.
    files : list
        List of paths with the PNG files.
    fps : int (optional)
        Frames per second.
    loop : int
        The number of times the GIF should loop.
        0 means that it will loop forever.
    """
    imgs = [Image.open(file) for file in files]
    imgs[0].save(fp=outfile, format='GIF', append_images=imgs[1:],
                 save_all=True, duration=int(1000/fps), loop=loop)


def plot_panels(nodes, panels, ax=None, plot_nodes=False, **plot_kwargs):
    '''The function `plot_panels` plots origami panels in 3D as a 3D collection of polygons using
    Matplotlib.
    
    Parameters
    ----------
    nodes
        The `nodes` parameter is an ndarray containing the coordinates of the vertices in 3D space. Each
    row in the array represents the coordinates of a single vertex, and the shape of the array is
    `(n_nodes, 3)`, where `n_nodes` is the number of vertices.
    panels
        The `panels` parameter is a list that contains the vertices number for each panel. Each element in
    the list represents a panel, and the numbers within that element correspond to the vertices that
    make up that panel.
    ax
        The `ax` parameter in the `plot_panels` function is a Matplotlib.axes object. It is used to specify
    the axes where the 3D plot of origami panels will be added. If `ax` is not provided (i.e., it is
    `None`), the function
    
    Returns
    -------
        The function `plot_panels` returns the Matplotlib axes object `ax` where the 3D graphic of the
    origami panels is plotted.
    
    '''
    """Plot origami panels in 3D as a 3D collection of polygons
    
    Parameters
    ----------
    nodes : ndarray, float
        Coordinates of the vertices (n_nodes, 3).
    panels : list
        List with the vertices number for each panel.
        The numbers should be integers.
    ax : Matplotlib.axes (optional)
        Axes to add the graphic. None by default. If None is
        passed it creates a new one.
    fix_aspect_ratio : bool
        Flag to fix the aspect ratio of the figure according
        to the location of the nodes. True by default.
    Returns
    -------
    ax : Matplotlib.axesg
        Axes to add the graphic. If None is passed it creates a
        new one.
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
    poly3d = list(nodes[panel] for panel in panels)
    ls = LightSource()
    poly_collection = Poly3DCollection(poly3d, 
                                       shade=True,
                                       lightsource=ls,
                                       **plot_kwargs)
    ax.add_collection3d(poly_collection)


    if plot_nodes==True:
        # Extract coordinates
        X = [n[0] for n in nodes]
        Y = [n[1] for n in nodes]
        Z = [n[2] for n in nodes]
        
        # Plot the nodes
        ax.scatter(X, Y, Z, c='blue', marker='o')
        
        # Annotate each node with its index
        for i, (x, y, z) in enumerate(nodes):
            ax.text(x, y, z, '%d' % i, size=10, zorder=1, color='red')
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_box_aspect([1,1,2])
    
    
    return ax


if __name__ == "__main__":
    
    # Example with 3 panels
    
    nodes = np.array([
        [0.000, 0, 0],
        [0.707, 0, 0.707],
        [1.414, 0, 0],
        [0.000, 1, 0],
        [0.707, 1, 0.707],
        [1.414, 1, 0.707]])
    

    
    
    # One square and two triangles
    panels = [
        [0, 1, 4, 3],
        [1, 2, 4],
        [2, 5, 4]]
    


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    plot_kwargs_wire = {"linewidths": 1, "edgecolors": "#3c3c3c",
                   "facecolors": "#d86a96", "alpha": 0.4}
    plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
                   "facecolors": "#d86a96"}
    nodes2 = nodes.copy()
    niter = 20
    dz = 1/niter
    files = []
    for cont in range(niter):
        
        plt.cla()
        
        nodes2[:, 2] += dz
        
        plot_panels(nodes, panels, ax=ax, **plot_kwargs_wire)
        plot_panels(nodes2, panels, ax=ax, **plot_kwargs)
        ax.auto_scale_xyz([0, 1.5], [0, 1], [0, 2])
        plt.axis("image")
        file = f"ori_{str(cont).zfill(2)}.png"
        plt.savefig(file)
        files.append(file)

    plt.show()

    save_gif_PIL("ori_anim.gif", files, fps=5, loop=0)

    [os.remove(file) for file in files]







# # -*- coding: utf-8 -*-
# """
# Visualize origami panels in Matplotlib

# @author: Nicolás Guarín-Zapata
# @date: May 2024
# """
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# from matplotlib.colors import LightSource
# import time



# def plot_panels(nodes, panels, ax=None, fix_aspect_ratio=True):
#     """Plot origami panels in 3D as a 3D collection of polygons
    
#     Parameters
#     ----------
#     nodes : ndarray, float
#         Coordinates of the vertices (n_nodes, 3).
#     panels : list
#         List with the vertices number for each panel.
#         The numbers should be integers.
#     ax : Matplotlib.axes (optional)
#         Axes to add the graphic. None by default. If None is
#         passed it creates a new one.
#     fix_aspect_ratio : bool
#         Flag to fix the aspect ratio of the figure according
#         to the location of the nodes. True by default.

#     Returns
#     -------
#     ax : Matplotlib.axesg
#         Axes to add the graphic. If None is passed it creates a
#         new one.
#     """
#     if ax is None:
#         fig = plt.figure()
#         ax = fig.add_subplot(111, projection='3d')
        
        
#     if fix_aspect_ratio:
#         x, y, z = nodes.T
#         max_range = np.array([x.max()-x.min(), y.max()-y.min(),
#                               z.max()-z.min()]).max() / 2.0
#         mean_x = x.mean()
#         mean_y = y.mean()
#         mean_z = z.mean()
#         ax.set_xlim(mean_x - max_range, mean_x + max_range)
#         ax.set_ylim(mean_y - max_range, mean_y + max_range)
#         ax.set_zlim(mean_z - max_range, mean_z + max_range)
        
        
#     poly3d = list(nodes[panel] for panel in panels)
#     ls = LightSource()
#     poly_collection = Poly3DCollection(poly3d, linewidths=1,
#                                        facecolors='#d86a96', 
#                                        edgecolors="#3c3c3c",
#                                        shade=True,
#                                        lightsource=ls)
#     ax.add_collection3d(poly_collection)
    


    
#     return ax

    
# if __name__ == "__main__":
    
#     # Example with 3 panels
    
#     nodes = np.array([
#         [0.000, 0, 0],
#         [0.707, 0, 0.707],
#         [1.414, 0, 0],
#         [0.000, 1, 0],
#         [0.707, 1, 0.707],
#         [1.414, 1, 0.707]])
    
#     # One square and two triangles
#     panels = [
#         [0, 1, 4, 3],
#         [1, 2, 4],
#         [2, 5, 4]]
    
    
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     plot_panels(nodes, panels, ax=ax, fix_aspect_ratio=True)
    
#     plt.show()