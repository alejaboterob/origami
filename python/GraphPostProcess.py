import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, writers
# from PIL import Image
# from PlotOri import PlotOri

def GraphPostProcess(U_his, STAT):
    '''The function `GraphPostProcess` generates a plot showing different types of stored energy over
    pseudo-time increments.
    
    Parameters
    ----------
    U_his
        It seems like you were about to provide some information about the `U_his` parameter, but the
    message got cut off. Could you please provide more details or let me know how I can assist you
    further with this code snippet?
    STAT
        STAT is a dictionary containing energy statistics for different components of a system. It includes
    keys such as 'PE' for total profile energy, 'bend' for bending energy, and 'bar' for stretching
    energy of bars. The values associated with these keys are arrays representing the energy values at
    each increment
    
    '''

    fig = plt.figure()

    plt.plot(np.arange(1, U_his.shape[1] + 1), STAT['PE']['strain'], color='#028090', linewidth=2)  # Total profile
    plt.plot(np.arange(1, U_his.shape[1] + 1), STAT['bend']['PE'] + STAT['bar']['PE'], color='#fe5d26')  # Folding energy
    plt.plot(np.arange(1, U_his.shape[1] + 1), STAT['bar']['PE'], color='#e574bc')  # Stretching energy of bars
    plt.xlabel('Increment Number (Pseudo-time)', fontsize=14)
    plt.ylabel('Stored Energy', fontsize=14)
    plt.grid(True)
    plt.show(block = False)

