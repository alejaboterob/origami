import numpy as np

def ConfigYoshimura(sec_hor, sec_vert, radius, height, fdang):
    """
    Configures the Yoshimura pattern with given parameters.
    
    sec_hor: Number of horizontal sections (divisions along the circumference)
    sec_vert: Number of vertical sections (divisions along the height)
    radius: Radius of the cylinder
    height: Height of the cylinder
    fdang: Folding angle (angle of deformation)

    Returns:
    Node: List of nodes (vertices)
    Panel: List of panels (faces)
    BDRY: Boundary edges or nodes
    """
    Node = []
    Panel = []
    BDRY = []

    # Angle and height step size
    d_theta = 2 * np.pi / sec_hor
    d_height = height / sec_vert

    X, Y, Z = [], [], []
    
    # Generate nodes (vertices)
    for i in range(sec_vert + 1):
        for j in range(sec_hor):
            x = radius * np.cos(j * d_theta)
            y = radius * np.sin(j * d_theta)
            z = i * d_height
            X.append(x)
            Y.append(y)
            Z.append(z)

    Node = np.column_stack((X, Y, Z))

    # Generate panels (faces)
    for i in range(sec_vert):
        for j in range(sec_hor):
            next_j = (j + 1) % sec_hor
            
            v1 = i * sec_hor + j
            v2 = i * sec_hor + next_j
            v3 = (i + 1) * sec_hor + j
            v4 = (i + 1) * sec_hor + next_j
            
            if (i + j) % 2 == 0:
                # Diagonal from v1 to v4
                Panel.append([v1, v2, v4])
                Panel.append([v1, v4, v3])
            else:
                # Diagonal from v2 to v3
                Panel.append([v1, v2, v3])
                Panel.append([v3, v4, v2])
    
    # Generate boundaries (optional)
    for j in range(sec_hor):
        BDRY.append(j)
        BDRY.append(sec_vert * sec_hor + j)

    return Node, Panel, BDRY


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def rotate_vertex_around_axis(vertex, axis_point1, axis_point2, angle):
    """
    Rotates a vertex around a specified axis defined by two points (axis_point1, axis_point2).
    
    vertex: The vertex to be rotated.
    axis_point1, axis_point2: Two points defining the rotation axis.
    angle: The angle to rotate the vertex by (fdang).
    
    Returns the rotated vertex coordinates.
    """
    # Translate vertex and axis points so axis_point1 is at the origin
    v = vertex - axis_point1
    p2 = axis_point2 - axis_point1
    
    # Normalize the axis direction
    axis = p2 / np.linalg.norm(p2)
    
    # Rodrigues' rotation formula
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    cross_prod = np.cross(axis, v)
    dot_prod = np.dot(axis, v)
    
    rotated_vertex = (v * cos_angle + cross_prod * sin_angle + axis * dot_prod * (1 - cos_angle))
    
    # Translate the vertex back to the original position
    return rotated_vertex + axis_point1

def ConfigYoshimura(sec_hor, sec_vert, radius, height, fold_fraction):
    """
    Configures the Yoshimura pattern with given parameters and applies folding.
    
    sec_hor: Number of horizontal sections (divisions along the circumference)
    sec_vert: Number of vertical sections (divisions along the height)
    radius: Radius of the cylinder
    height: Height of the cylinder
    fold_fraction: Fraction of the full folding (0 to 1, where 0 is unfolded and 1 is fully folded)
    
    Returns:
    Node: List of nodes (vertices)
    Panel: List of panels (faces)
    """
    Node = []
    Panel = []
    BDRY = []

    # Angle and height step size
    d_theta = 2 * np.pi / sec_hor
    d_height = height / sec_vert
    
    # Folding angle
    max_fdang = np.pi / 2  # Maximum folding angle for full fold
    fdang = max_fdang * fold_fraction
    
    # Generate nodes (vertices)
    for i in range(sec_vert + 1):
        for j in range(sec_hor):
            x = radius * np.cos(j * d_theta)
            y = radius * np.sin(j * d_theta)
            z = i * d_height
            Node.append(np.array([x, y, z]))
    
    # Generate panels (faces) and apply folding
    for i in range(sec_vert):
        for j in range(sec_hor):
            next_j = (j + 1) % sec_hor
            
            v1 = i * sec_hor + j
            v2 = i * sec_hor + next_j
            v3 = (i + 1) * sec_hor + j
            v4 = (i + 1) * sec_hor + next_j
            
            # First triangle (v1, v2, v4)
            Node[v4] = rotate_vertex_around_axis(Node[v4], Node[v1], Node[v2], fdang)
            Panel.append([v1, v2, v4])
            
            # Second triangle (v1, v4, v3)
            Node[v3] = rotate_vertex_around_axis(Node[v3], Node[v1], Node[v4], fdang)
            Panel.append([v1, v4, v3])
    
    return Node, Panel

def plot_yoshimura(Node, Panel):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for panel in Panel:
        poly3d = [[Node[vertex] for vertex in panel]]
        ax.add_collection3d(Poly3DCollection(poly3d, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25))
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1,1,2])  # Aspect ratio is 1:1:2
    
    plt.show()

# Example usage
sec_hor = 4  # Number of sections around the circumference
sec_vert = 2  # Number of sections along the height
radius = 1  # Radius of the cylinder
height = 2  # Height of the cylinder

# Simulate folding progression from unfolded to fully folded
for fold_fraction in [0.0, 0.25, 0.5, 0.75, 1.0]:
    Node, Panel = ConfigYoshimura(sec_hor, sec_vert, radius, height, fold_fraction)
    plot_yoshimura(Node, Panel)

