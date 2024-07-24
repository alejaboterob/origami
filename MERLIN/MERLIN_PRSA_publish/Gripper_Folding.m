%% =========== VALLEY FOLDING =========================================== %%
clear all; close all; clc;
%% Define geometry and material 

a = 4;
b = 2;
c = 6;
d = 2;
e = 2;

Node = [
    0, -b/2, 0;
    a, -b/2, 0;
    0, b/2, 0;
    a, b/2, 0;
    (a-e)/2, -(d+b)/2, -d/2;
    (a+e)/2, -(d+b)/2, -d/2;
    (a-e)/2, (d+b)/2, -d/2;
    (a+e)/2, (d+b)/2, -d/2;
    0, -b/2, c;
    0, b/2, c;
    a, -b/2, c;
    a, b/2, c
];

Panel = {
    [0, 1, 3, 2]+1;
    [0, 1, 5, 4]+1;
    [2, 3, 7, 6]+1;
    [0, 4, 8]+1;
    [0, 2, 9, 8]+1;
    [2, 6, 9]+1;
    [1, 10, 5]+1;
    [1, 3, 11, 10]+1;
    [3, 7, 11]+1
};

BDRY = [
    8, 4;
    4, 5;
    5, 10;
    10, 11;
    11, 7;
    7, 6;
    6, 9;
    9, 8
]+1;


%%
MaxIcr = 60; blam = 0.01; 
% Maximum increment number & initial load factor
Kf = 1e-1; Kb = Kf*1e5; E0 = 1e6; Abar = 1e-1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 0.1; limrht = 180-0.1;
% Left and right limits for the linear range of rotational stiffness
%%
BarMater = @(Ex)Ogden(Ex, E0); % Define bar material constitutive
RotSpring = @(he,h0,kpi,L0)EnhancedLinear(he,h0,kpi,L0,limlft,limrht);

%% Set up boundary conditions
% Set up boundary conditions
Supp = [
    1, 1, 1, 1;
    2, 1, 1, 1;
    3, 1, 1, 1;
    4, 1, 1, 1
];

indp = [5, 6, 7,8];
ff = [1, 1, -1, -1] .*3;  
Load = [indp', zeros(length(indp), 1), ff', zeros(length(indp), 1)];
indp = Load(:, 1);

%% Perform analysis
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
truss.U0 = zeros(3*size(truss.Node,1),1);
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr);
U_his = real(U_his);
LF_his = real(LF_his);

%% Visualize simulation
instdof = -(indp(1)*3-2);
interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','miura5x5fold',0.05,LF_his,instdof,[-inf inf -inf inf])
% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','miura5x5fold',0.0001)

%% Plot stored energy vs. pseudo time
% Red line is the total profile. Between red and cyan is the folding
% energy. Between cyan and magenta is the portion of energy for bending. 
% Below magenta is the stretching energy of bars.
STAT = PostProcess(Data,truss,angles);  
figure()
plot(1:size(U_his,2),STAT.PE.strain,'r-','linewidth',2);
grid on
hold on
plot(1:size(U_his,2),STAT.bend.PE+STAT.bar.PE,'c-');
plot(1:size(U_his,2),STAT.bar.PE,'m-');
xlabel('Icrement Number (Pseudo-time)','fontsize',14);
ylabel('Stored Energy','fontsize',14);

%% Plot final configuration
Ux = U_his(:,end);
Nodew = truss.Node;
Nodew(:,1) = truss.Node(:,1)+Ux(1:3:end);
Nodew(:,2) = truss.Node(:,2)+Ux(2:3:end);
Nodew(:,3) = truss.Node(:,3)+Ux(3:3:end); 
figure()
PlotOri(Nodew,angles.Panel,truss.Trigl);
axis equal; axis off;
camproj('perspective')
light
view(117,18)
rotate3d on