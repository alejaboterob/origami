# -*- coding: utf-8 -*-
"""
Visualize origami panels in Matplotlib

@author: Nicolás Guarín-Zapata
@date: May 2024
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import LightSource
import time



def plot_panels(nodes, panels, ax=None, fix_aspect_ratio=True):
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
        
        
    if fix_aspect_ratio:
        x, y, z = nodes.T
        max_range = np.array([x.max()-x.min(), y.max()-y.min(),
                              z.max()-z.min()]).max() / 2.0
        mean_x = x.mean()
        mean_y = y.mean()
        mean_z = z.mean()
        ax.set_xlim(mean_x - max_range, mean_x + max_range)
        ax.set_ylim(mean_y - max_range, mean_y + max_range)
        ax.set_zlim(mean_z - max_range, mean_z + max_range)
        
        
    poly3d = list(nodes[panel] for panel in panels)
    ls = LightSource()
    poly_collection = Poly3DCollection(poly3d, linewidths=1,
                                       facecolors='#d86a96', 
                                       edgecolors="#3c3c3c",
                                       shade=True,
                                       lightsource=ls)
    ax.add_collection3d(poly_collection)
    


    
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
    plot_panels(nodes, panels, ax=ax, fix_aspect_ratio=True)
    
    plt.show()