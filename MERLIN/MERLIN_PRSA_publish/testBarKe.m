Ui_example = [0.11449873, 0.60748948, 0.92399539, 0.02126101, 0.24440479, 0.48081165]'
B_example = [0.99377798, 0.51655026, 0.52905187, 0.89969584, 0.50178366, 0.91173732]
L_example = 10  
E0 = 1e6
CM_example = @(Ex)Ogden(Ex, E0); % Define bar material constitutive
A_example = 1.0  
[~,Rbe,Kbe] = BarKe(Ui_example,B_example,L_example,CM_example,A_example)
