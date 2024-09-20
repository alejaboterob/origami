import json
import numpy as np


def readFOLD(file_path):
    try:
        # Open and read the file as a regular text file
        with open(file_path, 'r') as file:
            file_content = file.read()  # Read the entire content of the file

        # Attempt to parse the content as JSON
        data_dict = json.loads(file_content)  # This loads the JSON content into a dictionary
        print(data_dict)  # Output the dictionary to confirm it's loaded correctly

    except FileNotFoundError:
        print(f"File not found: {file_path}")

    Nodes = np.array(data_dict['vertices_coords'])
    Panels = np.array(data_dict['faces_vertices'], dtype=object)
    fda = np.array(data_dict['edges_foldAngle']) 
    BDRY = np.where(np.array(fda) == None)

    return Nodes, Panels, BDRY