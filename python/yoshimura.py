import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_yoshimura_pattern(radius, height, n_circumference, n_height):
    """
    Generates the vertices and panels for a Yoshimura pattern on a cylinder.
    
    radius: Radius of the cylinder
    height: Height of the cylinder
    n_circumference: Number of divisions along the circumference
    n_height: Number of divisions along the height
    
    Returns vertices and panels.
    """
    vertices = []
    panels = []
    
    # Angle and height step size
    d_theta = 2 * np.pi / n_circumference
    d_height = height / n_height
    
    # Generate vertices
    for i in range(n_height + 1):
        for j in range(n_circumference):
            x = radius * np.cos(j * d_theta)
            y = radius * np.sin(j * d_theta)
            z = i * d_height
            vertices.append([x, y, z])
    
    # Generate panels
    for i in range(n_height):
        for j in range(n_circumference):
            next_j = (j + 1) % n_circumference
            
            v1 = i * n_circumference + j
            v2 = i * n_circumference + next_j
            v3 = (i + 1) * n_circumference + j
            v4 = (i + 1) * n_circumference + next_j
            
            if (i + j) % 2 == 0:
                # Diagonal from v1 to v4
                panels.append([vertices[v1], vertices[v2], vertices[v4]])
                panels.append([vertices[v1], vertices[v4], vertices[v3]])
            else:
                # Diagonal from v2 to v3
                panels.append([vertices[v1], vertices[v2], vertices[v3]])
                panels.append([vertices[v3], vertices[v4], vertices[v2]])
    
    return vertices, panels

# Define parameters
radius = 1  # Radius of the cylinder
height = 2  # Height of the cylinder
n_circumference = 16  # Number of divisions along the circumference
n_height = 8  # Number of divisions along the height

# Generate pattern
vertices, panels = generate_yoshimura_pattern(radius, height, n_circumference, n_height)


def plot_yoshimura_pattern(panels):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for panel in panels:
        x = [vertex[0] for vertex in panel]
        y = [vertex[1] for vertex in panel]
        z = [vertex[2] for vertex in panel]
        
        # Close the loop by adding the first point to the end
        x.append(panel[0][0])
        y.append(panel[0][1])
        z.append(panel[0][2])
        
        ax.plot(x, y, z, 'bo-')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    plt.show()

# Plot the Yoshimura pattern
plot_yoshimura_pattern(panels)
