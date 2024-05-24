% Load variables from the .mat file saved by Python
python_variables = load('C:\Users\mboter45\OneDrive - Universidad EAFIT\Estructuras origami\proyecto\code\Paulino\python\variables.mat');

% % Compare variables with MATLAB variables
% Assuming you already have the variables NODE, PANEL, and BDRY in MATLAB
Panel = cell2mat(Panel);

isequal(python_variables.NODE, Node)
isequal(python_variables.PANEL, Panel)
isequal(python_variables.BDRY, BDRY)

% isequal(NODE2, Node)
% isequal(PANEL2, Panel)
% isequal(BDRY2, BDRY)
