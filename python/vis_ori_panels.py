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


nodes = np.array([
    [0.000, 0, 0],
    [0.707, 0, 0.707],
    [1.414, 0, 0],
    [0.000, 1, 0],
    [0.707, 1, 0.707],
    [1.414, 1, 0]])

panels = [
    [0, 1, 4, 3],
    [1, 2, 5, 4]]

x, y, z = nodes.T

poly3d = [nodes[panels[0]], nodes[panels[1]]]

#%%
ls = LightSource()
poly_collection = Poly3DCollection(poly3d, linewidths=1,
                                   facecolors='#d86a96', 
                                   edgecolors="#3c3c3c", shade=True,
                                   lightsource=ls)


#%% Visualization
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.add_collection3d(poly_collection)


# Fix aspect ratio
max_range = np.array([x.max()-x.min(), y.max()-y.min(),
                      z.max()-z.min()]).max() / 2.0
mean_x = x.mean()
mean_y = y.mean()
mean_z = z.mean()
ax.set_xlim(mean_x - max_range, mean_x + max_range)
ax.set_ylim(mean_y - max_range, mean_y + max_range)
ax.set_zlim(mean_z - max_range, mean_z + max_range)

plt.show()