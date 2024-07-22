import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, writers
# from PIL import Image
# from PlotOri import PlotOri
from plot_ori_panels import plot_panels, save_gif_PIL

files = []

def displacement(U_his, truss, angles, instdof, pausetime, LF_his=None):
    '''The `displacement` function visualizes the displacement and load factor history of a truss structure
    over time.

    Parameters
    ----------
    U_his
        The `U_his` parameter seems to represent the history of displacements at specific degrees of
    freedom in a structural analysis context. It likely contains the displacement values over time for
    different degrees of freedom of a structure under loading conditions.
    truss
        The `truss` parameter in the `displacement` function is likely a data structure representing a
    truss system. It could contain information such as node coordinates, element connectivity, material
    properties, and boundary conditions for the truss structure. This information is essential for
    calculating the displacement and visualizing the
    angles
        The `angles` parameter in the `displacement` function seems to be unused in the provided code
    snippet. If you need assistance with how to use this parameter or have any specific questions
    related to it, please let me know.
    instdof
        The `instdof` parameter in the `displacement` function seems to represent the degree of freedom for
    which the displacement is being calculated. It is used to extract the displacement values from the
    `U_his` array based on the index provided by `instdof`.
    pausetime
        The `pausetime` parameter in the `displacement` function is used to specify the time interval (in
    seconds) to pause between each frame of the animation. This allows you to control the speed at which
    the animation is displayed to the user. You can adjust this parameter based on how fast
    LF_his
        The `LF_his` parameter in the `displacement` function seems to represent the history of load
    factors. It is used in the function to plot the load factor values corresponding to the
    displacements over time. The function iterates through the load factor history and updates the plot
    with each iteration.

    '''

    fig = plt.figure(figsize=(8, 6))
    fig.set_facecolor('w')
    dsp = np.sign(instdof) * U_his[np.abs(instdof).astype('int')]

    for i in range(len(LF_his)):
        plt.cla()
        plt.plot(dsp[:i+1], LF_his[:i+1], color='#e574bc', linewidth=2)
        plt.plot(dsp[i], LF_his[i], 'o', color='#028090', linewidth=2)
        plt.xlabel('displacement', fontsize=14)
        plt.ylabel('load factor', fontsize=14)
        plt.pause(pausetime)


    plt.show(block = False)
    # save_gif_PIL("ori_anim.gif", files, fps=5, loop=0)
