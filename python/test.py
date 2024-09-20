import numpy as np
import json

def construct_waterbomb_pattern(rows, cols, cell_width, cell_height):
    """
    Constructs a waterbomb pattern with specified rows and columns.

    Parameters:
    - rows: Number of rows in the pattern.
    - cols: Number of columns in the pattern.
    - cell_width: Width of each cell in the pattern.
    - cell_height: Height of each cell in the pattern.

    Returns:
    - vertices: List of vertex coordinates.
    - edges: List of edges connecting vertices.
    - edge_types: List of edge types ('mountain' or 'valley').
    """
    vertices = []
    edges = []
    edge_types = []

    # Generate vertices
    for i in range(rows + 1):
        for j in range(cols + 1):
            x = j * cell_width
            y = i * cell_height
            vertices.append([x, y, 0.0])  # Assuming z=0 for 2D pattern

    # Generate edges
    for i in range(rows):
        for j in range(cols):
            # Define vertex indices for each square cell
            top_left = i * (cols + 1) + j
            top_right = top_left + 1
            bottom_left = (i + 1) * (cols + 1) + j
            bottom_right = bottom_left + 1

            # Add horizontal edges
            edges.append([top_left, top_right])
            edge_types.append('valley' if (i + j) % 2 == 0 else 'mountain')

            edges.append([bottom_left, bottom_right])
            edge_types.append('valley' if (i + j) % 2 == 0 else 'mountain')

            # Add vertical edges
            edges.append([top_left, bottom_left])
            edge_types.append('valley' if (i + j) % 2 == 0 else 'mountain')

            edges.append([top_right, bottom_right])
            edge_types.append('valley' if (i + j) % 2 == 0 else 'mountain')

            # Add diagonal edges for the waterbomb folds
            edges.append([top_left, bottom_right])
            edge_types.append('mountain' if (i + j) % 2 == 0 else 'valley')

            edges.append([top_right, bottom_left])
            edge_types.append('mountain' if (i + j) % 2 == 0 else 'valley')

    return vertices, edges, edge_types

# Example usage
rows = 4  # Number of rows
cols = 4  # Number of columns
cell_width = 1.0  # Width of each cell
cell_height = 1.0  # Height of each cell

vertices, edges, edge_types = construct_waterbomb_pattern(rows, cols, cell_width, cell_height)

# Save pattern to JSON format
pattern_data = {
    'vertices_coords': vertices,
    'edges_vertices': edges,
    'edges_assignment': edge_types
}

