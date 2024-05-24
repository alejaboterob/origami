import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
from PIL import Image
from PlotOri import PlotOri

def VisualFold(U_his, truss, angles, recordtype, filename, pausetime, LF_his=None, instdof=None, axislim=None):
    """
    Record the simulation to file if needed:
    recordtype = 'none': do not save the simulation
    recordtype = 'video': save simulation in MP4 format;
    recordtype = 'imggif': save simulatrion in GIF format.
    pausement: pause time between each frame (in seconds).
               If recordtype = 'none': use a small number such as 0.0001;
               Otherwise, pausetime = 1/fps;
    If input data does not include 'load_his', 'instdof', 'axislim', the
    function does not plot load vs. displacement diagram.
    axislim: Axis limits (bounding box) for load vs. displacement diagram.
             format: [xmin,xmax,ymin,ymax].
    instdof: specify the DOF of interest for displacement measure.
    """
    Node = truss['Node']
    Trigl = truss['Trigl']
    Panel = angles['Panel']

    if LF_his is not None and len(LF_his.shape) > 1:
        LF_his = np.sum(LF_his, axis=1)

    if recordtype == 'video':
        plt.rcParams['animation.ffmpeg_path'] = '/path/to/ffmpeg'
        writer = writers['ffmpeg'](fps=1/pausetime)
        fig = plt.figure('units', 'pixels')
        fig.set_facecolor('w')
    elif recordtype == 'imggif':
        filename = f"{filename}.gif"
        fig = plt.figure('units', 'pixels')
        fig.set_facecolor('w')
    else:
        print('Not recording')
        fig = plt.figure()
        # fig = plt.figure('units', 'pixels')
        # fig.set_facecolor('w')

    plt.ion()
    for i in range(U_his.shape[1]):
        U = U_his[:, i]
        plt.clf()
        # plt.view(35, 30)
        Nodew = Node.copy()
        Nodew[:, 0] = Node[:, 0] + U[::3]
        Nodew[:, 1] = Node[:, 1] + U[1::3]
        Nodew[:, 2] = Node[:, 2] + U[2::3]

        PlotOri(Node, Panel, Trigl, '-', 0.3, 'none')
        PlotOri(Nodew, Panel, Trigl)
        plt.axis('equal')
        plt.axis('off')
        plt.light()
        plt.pause(pausetime)

        if recordtype == 'imggif':
            frame = cv2.VideoCapture(f2).read()[1]
            im = Image.fromarray(frame)
            imind, cm = im.convert('P', palette=Image.ADAPTIVE, colors=256).getpalette()
            if i == 1:
                im.save(filename, 'GIF', save_all=True, append_images=[], duration=0, loop=0)
            else:
                im.save(filename, 'GIF', save_all=True, append_images=[], duration=pausetime, loop=0)
        elif recordtype == 'video':
            frame = cv2.VideoCapture(f2).read()[1]
            writerObj.write(frame)

            cv2.waitKey(0)
            cv2.destroyAllWindows()
        if recordtype == 'video':
            writerObj.release()

    plt.ioff()
    plt.close(fig)

    if recordtype == 'video':
        writer.finish()
