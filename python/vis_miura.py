# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:28:44 2024

@author: nguarinz
"""
import numpy as np
import matplotlib.pyplot as plt
from plot_ori_panels import plot_panels, save_gif_PIL

data = np.load("data.npy", allow_pickle=True).item()
nodes_ini = data["Node"]
nodes_end = data["Nodew"]
panels = data["Panel"]

mins = np.minimum(nodes_ini.min(axis=0), nodes_end.min(axis=0))
maxs = np.maximum(nodes_ini.max(axis=0), nodes_end.max(axis=0))

plot_kwargs_wire = {"linewidths": 1, "edgecolors": "#3c3c3c22",
               "facecolors": "#d86a96", "alpha": 0.0, "zorder": 4,
               "zsort": "max"}
plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
               "facecolors": "#d86a96", "zorder": 5,
               "zsort": "min"}
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot_panels(nodes_ini, panels, ax=ax, **plot_kwargs_wire)
plot_panels(nodes_end, panels, ax=ax, **plot_kwargs)
plt.xticks([0, 5, 10, 15, 20])
plt.yticks([0, 5, 10, 15, 20])
ax.set_zticks([0, 2])
ax.set(xticklabels=[], yticklabels=[], zticklabels=[])
ax.auto_scale_xyz([mins[0], maxs[0]],
                  [mins[1], maxs[1]],
                  [mins[2], maxs[2]])
plt.axis("image")