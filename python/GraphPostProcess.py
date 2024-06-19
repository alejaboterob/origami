import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, writers
# from PIL import Image
# from PlotOri import PlotOri

def GraphPostProcess(U_his, STAT):

    fig = plt.figure()

    plt.plot(np.arange(1, U_his.shape[1] + 1), STAT['PE']['strain'], color='#028090', linewidth=2)  # Total profile
    plt.plot(np.arange(1, U_his.shape[1] + 1), STAT['bend']['PE'] + STAT['bar']['PE'], color='#fe5d26')  # Folding energy
    plt.plot(np.arange(1, U_his.shape[1] + 1), STAT['bar']['PE'], color='#e574bc')  # Stretching energy of bars
    plt.xlabel('Increment Number (Pseudo-time)', fontsize=14)
    plt.ylabel('Stored Energy', fontsize=14)
    plt.grid(True)
    plt.show()

