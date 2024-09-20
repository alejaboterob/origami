import json
import numpy as np

def configure_waterbomb_pattern(data, num_horizontal_sections, num_vertical_sections, length_horizontal, length_vertical):
    """
    Configure a waterbomb pattern with specified sections and lengths.
    
    Parameters:
    - data: Original JSON data of the waterbomb pattern.
    - num_horizontal_sections: Number of sections along the horizontal axis.
    - num_vertical_sections: Number of sections along the vertical axis.
    - length_horizontal: Length of each section along the horizontal axis.
    - length_vertical: Length of each section along the vertical axis.
    
    Returns:
    - Modified JSON data with the new configuration.
    """
    
    # Calculate the new vertices coordinates
    new_vertices = []
    for i in range(num_vertical_sections + 1):
        for j in range(num_horizontal_sections + 1):
            x = j * length_horizontal
            y = i * length_vertical
            new_vertices.append([x, y, 0.0])  # Z coordinate is 0 for a flat waterbomb base

    # Create new edges based on the new grid
    new_edges = []
    for i in range(num_vertical_sections):
        for j in range(num_horizontal_sections):
            # Define vertices for each square in the grid
            top_left = i * (num_horizontal_sections + 1) + j
            top_right = top_left + 1
            bottom_left = (i + 1) * (num_horizontal_sections + 1) + j
            bottom_right = bottom_left + 1
            
            # Add horizontal and vertical edges
            new_edges.append([top_left, top_right])
            new_edges.append([top_left, bottom_left])
            new_edges.append([bottom_left, bottom_right])
            new_edges.append([top_right, bottom_right])

    # Update the JSON data with new vertices and edges
    data['vertices_coords'] = new_vertices
    data['edges_vertices'] = new_edges
    
    return data

# Example usage
# Load the original JSON data
file_path = r"C:\Users\mboter45\Downloads\waterbomb_10PercentFolded.fold"
with open(file_path, 'r') as file:
    data = json.load(file)

# Configure with desired parameters
num_horizontal_sections = 4
num_vertical_sections = 4
length_horizontal = 1.0
length_vertical = 1.0

modified_data = configure_waterbomb_pattern(data, num_horizontal_sections, num_vertical_sections, length_horizontal, length_vertical)











import json
import numpy as np

from plot_ori_panels import plot_panels, save_gif_PIL
import matplotlib.pyplot as plt

from ConfigYoshimura import plot_yoshimura_nodes

def ConfigWaterbomb(file_path):
    try:
        # Open and read the file as a regular text file
        with open(file_path, 'r') as file:
            file_content = file.read()  # Read the entire content of the file

        # Attempt to parse the content as JSON
        try:
            data_dict = json.loads(file_content)  # This loads the JSON content into a dictionary
            print(data_dict)  # Output the dictionary to confirm it's loaded correctly
        except json.JSONDecodeError:
            print("Error decoding JSON. The file content might not be in JSON format.")

    except FileNotFoundError:
        print(f"File not found: {file_path}")

    Nodes = np.array(data_dict['vertices_coords'])
    Panels = np.array(data_dict['faces_vertices'], dtype=object)
    fda = np.array(data_dict['edges_foldAngle']) 
    BDRY = np.where(np.array(fda) == None)

    return Nodes, Panels, BDRY




# file_path = r"C:\Users\mboter45\Downloads\waterbomb_10PercentFolded.fold"

# Nodes, Panels, BDRY = ConfigWaterbomb(file_path)
# print(Nodes)
# print(Panels)


# try:
#     # Open and read the file as a regular text file
#     with open(file_path, 'r') as file:
#         file_content = file.read()  # Read the entire content of the file

#     # Print the content of the file
#     print(file_content)  # Output the content to confirm it's loaded correctly

# except FileNotFoundError:
#     print(f"File not found: {file_path}")

# plot_yoshimura_nodes(Nodes)

# plot_kwargs = {"linewidths": 1, "edgecolors": "#3c3c3c",
#                 "facecolors": "#d86a96"}

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# plot_panels(Nodes, Panels, ax, **plot_kwargs)

# mins = np.minimum(Nodes.min(axis=0), Nodes.min(axis=0))
# maxs = np.maximum(Nodes.max(axis=0), Nodes.max(axis=0))

# # plt.xticks([])
# # plt.yticks([])
# ax.auto_scale_xyz([mins[0], maxs[0]],
#             [mins[1], maxs[1]],
#             [mins[2], maxs[2]])    
# # ax.axis('equal')
# # ax.set(xticklabels=[], yticklabels=[], zticklabels=[])


# plt.show(block = False)

# plt.show()

# # import json

# # # Correct file path where your file is located
# # file_path = 'C:/Users/mboter45/OneDrive - Universidad EAFIT/Estructuras origami/fold/waterbomb_10PercentFolded.fold'

# # try:
# #     # Open and read the file as a regular text file
# #     with open(file_path, 'r') as file:
# #         file_content = file.read()  # Read the entire content of the file

# #     # Attempt to parse the content as JSON
# #     try:
# #         data_dict = json.loads(file_content)  # This loads the JSON content into a dictionary
# #         print(data_dict)  # Output the dictionary to confirm it's loaded correctly
# #     except json.JSONDecodeError:
# #         print("Error decoding JSON. The file content might not be in JSON format.")

# # except FileNotFoundError:
# #     print(f"File not found: {file_path}")
