import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from plot_ori_panels import plot_panels, save_gif_PIL

def ConfigYoshimura(sec_hor, sec_vert, radius, height, fdang):
    """
    Configures the Yoshimura pattern with given parameters and applies folding.
    
    sec_hor: Number of horizontal sections (divisions along the circumference)
    sec_vert: Number of vertical sections (divisions along the height)
    radius: Radius of the cylinder
    height: Height of the pattern
    fold_fraction: Fraction of the full folding (0 to 1, where 0 is unfolded and 1 is fully folded)
    
    Returns:
    Node: List of nodes (vertices)
    Panel: List of panels (faces)
    """
    Node = []
    Panel = []

    # Angle and height step size
    d_theta = 2 * np.pi / sec_hor
    
    # Generate nodes (vertices)
    for i in range(sec_vert + 1):
        for j in range(sec_hor):
            if i % 2 == 0:
                x = radius * np.cos(j * d_theta)
                y = radius * np.sin(j * d_theta)
                z = i * height
                Node.append(np.array([x, y, z]))
            else:
                x = radius * np.cos(j * d_theta + d_theta/2)
                y = radius * np.sin(j * d_theta + d_theta/2)
                z = i * height
                Node.append(np.array([x, y, z]))
    
    Node_idx = np.arange(len(Node)).reshape(sec_vert + 1, sec_hor)
    Panel = np.empty((0, 3))
    for i in range(sec_vert):
        if i % 2 == 0:
            Panel_tmp = np.vstack((
                np.vstack((Node_idx[i], 
                        np.concatenate((Node_idx[i, 1:], [Node_idx[i, 0]])), 
                        Node_idx[i + 1])).T ,
                np.vstack((Node_idx[i+1], 
                        np.concatenate((Node_idx[i+1,1:],[Node_idx[i+1,0]])), 
                        np.concatenate((Node_idx[i,1:],[Node_idx[i,0]])))).T ))
            Panel = np.append(Panel, Panel_tmp, axis=0)   
        else:
            Panel_tmp = np.vstack(( 
                np.vstack((Node_idx[i], 
                        np.concatenate((Node_idx[i, 1:], [Node_idx[i+1, 0]])), 
                        np.concatenate((Node_idx[i+1, 1:], [Node_idx[i, 0]])))).T,
                np.vstack((Node_idx[i], 
                        Node_idx[i+1], 
                        np.concatenate((Node_idx[i+1, 1:], [Node_idx[i+1, 0]])))).T )) 
            Panel = np.append(Panel, Panel_tmp, axis=0) 

    Node = np.array(Node)
    Panel = np.array(Panel, dtype=int)  

    
    return Node, Panel, Node_idx
    
    # Generate panels (faces)
    # for i in range(sec_vert):
    #     for j in range(sec_hor):
    #         idx = i * (sec_hor + 1) + j
    #         Panel.append([idx, idx + 1, idx + sec_hor])
             
    #         Panel.append([idx + 1, idx + sec_hor + 2, idx + sec_hor + 1])

    #         idx = (i - 1) * (sec_hor + 1) + j
    #         Panel.append([idx + sec_hor + 1, idx + sec_hor + 2, idx + 2 * sec_hor + 2])
    #         Panel.append([idx + sec_hor + 2, idx + 2 * sec_hor + 3, idx + 2 * sec_hor + 2])
    

def plot_yoshimura_nodes(Node):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Extract coordinates
    X = [n[0] for n in Node]
    Y = [n[1] for n in Node]
    Z = [n[2] for n in Node]
    
    # Plot the nodes
    ax.scatter(X, Y, Z, c='blue', marker='o')
    
    # Annotate each node with its index
    for i, (x, y, z) in enumerate(Node):
        ax.text(x, y, z, '%d' % i, size=10, zorder=1, color='red')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1,1,2])  # Aspect ratio is 1:1:2
    
    plt.show()


# Example usage
sec_hor = 4  # Number of sections around the circumference (minimum 3)
sec_vert = 8  # Number of sections along the height
radius = 1  # Radius of the cylinder
height = 1  # Height of the pattern

Node, Panel, Node_idx = ConfigYoshimura(sec_hor, sec_vert, radius, height, np.pi/8)
plot_yoshimura_nodes(Node)

Node = np.array(Node)
Panel = np.array(Panel, dtype=int)

plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
                "facecolors": "#d86a96"}

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot_panels(Node, Panel, ax, **plot_kwargs)

mins = np.minimum(Node.min(axis=0), Node.min(axis=0))
maxs = np.maximum(Node.max(axis=0), Node.max(axis=0))

# plt.xticks([])
# plt.yticks([])
ax.auto_scale_xyz([mins[0], maxs[0]],
            [mins[1], maxs[1]],
            [mins[2], maxs[2]])    
# ax.axis('equal')
# ax.set(xticklabels=[], yticklabels=[], zticklabels=[])


plt.show(block = False)

plt.show()
# save_gif_PIL(Node, Panel, 'yoshimura.gif', plot_kwargs, plot_panels, num_frames=100, fps=10)
