# -*- coding: utf-8 -*-
"""
Visualize origami triangular panels in Matplotlib

@author: Nicolás Guarín-Zapata
@date: May 2024
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource



nodes = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, -0.5],
    [0, 1, 0]])

panels = [
    [0, 1, 3],
    [1, 2, 3]]

x, y, z = nodes.T

fig = plt.figure()
ls = LightSource()
ax = fig.add_subplot(projection='3d')
ax.plot_trisurf(x, y, z, triangles=panels, linewidth=1,
                facecolor='#d86a96', edgecolor="#3c3c3c",
                shade=True, lightsource=ls)

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