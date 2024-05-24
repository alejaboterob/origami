import numpy as np
import matplotlib.pyplot as plt

def PlotOri(Node, Panel, Trigl, lstyle='-', alfa=1, color='g'):
    if len(Node) < 4:
        lstl = '-'
        al = 1
        cc = 'g'
    else:
        lstl = lstyle
        al = alfa
        cc = color
    
    # if Trigl.any():
        # plt.patch(faces=Trigl, vertices=Node, facecolor=cc, linestyle='none', facelighting='flat', edgecolor=(1-al)*np.ones(3))
    
    # plt.hold(True)
    
    Panelsize = [len(p) for p in Panel]
    Ptri = [Panel[i] for i in np.where(np.array(Panelsize)==3)[0]]
    Pquad = [Panel[i] for i in np.where(np.array(Panelsize)==4)[0]]
    
    # if Trigl.any():
        # plt.patch(faces=np.concatenate(Ptri), vertices=Node, facecolor='none', linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
        # plt.patch(faces=np.concatenate(Pquad), vertices=Node, facecolor='none', linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
    # else:
        # plt.patch(faces=np.concatenate(Ptri), vertices=Node, facecolor=cc, linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
        # plt.patch(faces=np.concatenate(Pquad), vertices=Node, facecolor=cc, linestyle=lstl, linewidth=1, edgecolor=(1-al)*np.ones(3))
    
    # plt.hold(False)

