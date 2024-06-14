import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, writers
# from PIL import Image
# from PlotOri import PlotOri
from plot_ori_panels import plot_panels, save_gif_PIL

files = []

def displacement(U_his, truss, angles, instdof, pausetime, LF_his=None):

    fig = plt.figure(figsize=(8, 6))
    fig.set_facecolor('w')
    dsp = np.sign(instdof) * U_his[np.abs(instdof).astype('int')]

    for i in range(len(LF_his)):
        plt.cla()
        plt.plot(dsp[:i+1], LF_his[:i+1], 'b-', linewidth=2)
        plt.plot(dsp[i], LF_his[i], 'ro', linewidth=2)
        plt.xlabel('displacement', fontsize=14)
        plt.ylabel('load factor', fontsize=14)
        plt.pause(pausetime)


    plt.show()
    # save_gif_PIL("ori_anim.gif", files, fps=5, loop=0)
