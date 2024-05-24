%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                  MERLIN                               %%
%                         Ke Liu, Glaucio H. Paulino                      %
% Ref: K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid    %
%      origami - An efficient computational approach.' Proceedings of     %
%      the Royal Society A.                                               %
%      K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to  % 
%      capture highly nonlinear behavior of non-rigid origami.'           %
%      Proceedings of IASS Annual Symposium 2016.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =========== MIURA BENDING =========================================== %%
clear all; close all; clc;
%% Define geomtry and material
sec_hor=5;  sec_vert=5;  % Number of unit cells in each direction
theta = 60; a = 2; b = 2; fdang = 54.7356; 
% Geometry of the Miura: a, b are the edge lengths of each parallelogram
% panel; fdang controls the folding angle; theta is the panel angle
MaxIcr = 98; blam = 0.03; 
% Maximum increment number & initial load factor
Kf = 1e-1; Kb = Kf*10; C0 = 1e4; Abar = 1e-1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 0.1; limrht = 360-0.1;
% Left and right limits for the linear range of rotational stiffness
[Node,Panel,~]=ConfigMiura(sec_hor,sec_vert,theta,a,b,fdang);
BarMater = @(Ex)Ogden(Ex, C0); % Define bar material constitutive
RotSpring = @(he,h0,kpi,L0)EnhancedLinear(he,h0,kpi,L0,limlft,limrht);
% Define rotational spring constitutive

%% Set up boundary conditions
m = size(Node,1);
Supp = [    56, 1, 0, 0;
            57, 1, 1, 1;
            65, 1, 0, 1];
indp = [56,66]';
ff = -1*ones(length(indp),1); 
Load = [indp,zeros(length(indp),1),zeros(length(indp),1),ff];
indp = Load(:,1);

%% Perform analysis
% Assemble input data
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
% Specify initial deformation state
truss.U0 = zeros(3*size(truss.Node,1),1);
% Perform path-following analysis using MGDCM
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr);
% Clean output data
U_his = real(U_his);
LF_his = real(LF_his);

%% Visualize simulation
instdof = -indp(1)*3;
interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','miura5x5bend',0.0001,LF_his,instdof,[-inf inf -inf inf])
% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','miura5x5bend',0.0001)

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
xlabel('Increment Number (Pseudo-time)','fontsize',14);
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