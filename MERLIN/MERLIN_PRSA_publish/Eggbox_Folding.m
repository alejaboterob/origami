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
%% =========== EGGBOX FOLDING =========================================== %%
clear all; close all; clc;
%% Define geomtry and material
sec_x=5;  sec_y=5;   % Number of unit cells in each direction
theta = 60; a = 1; b = 1; fdang = 45; 
% Geometry of the Eggbox: a, b are the edge lengths of each parallelogram
% panel; fdang controls the folding angle (<theta); theta is the panel angle
MaxIcr = 90; blam = 0.05; 
% Maximum increment number & initial load factor
Kf = 1e-1; Kb = Kf*1e4; E0 = 1e6; Abar = 1e-1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 0.1; limrht = 360-0.1;
% Left and right limits for the linear range of rotational stiffness
[Node,Panel,~]=ConfigEggbox(sec_x,sec_y,theta,a,b,fdang); %Eggbox
BarMater = @(Ex)Ogden(Ex, E0); % Define bar material constitutive
RotSpring = @(he,h0,kpi,L0)EnhancedLinear(he,h0,kpi,L0,limlft,limrht);
% Define rotational spring constitutive

%% Set up boundary conditions
leftx = 1:(2*sec_y+1);
leftz = 1:2:(2*sec_y+1);
rightz = [1:2:(2*sec_y+1)]+(2*sec_y+1)*(2*sec_x);
rightxp = [2:2:(2*sec_y+1)]+(2*sec_y+1)*(2*sec_x);
Supp = [6, 0, 1, 0;
        leftx',ones(numel(leftx),1),zeros(numel(leftx),1),zeros(numel(leftx),1);
        leftz',zeros(numel(leftz),1),zeros(numel(leftz),1),ones(numel(leftz),1);
%         rightz',zeros(numel(rightz),1),zeros(numel(rightz),1),ones(numel(rightz),1);
        rightxp(3), 0, 1, 1];
indp = rightxp';
ff = -1*ones(length(indp),1); 
Load = [indp,ff,zeros(length(indp),1),zeros(length(indp),1)];
indp = Load(:,1);

%% Perform analysis
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
truss.U0 = zeros(3*size(truss.Node,1),1);
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr);
U_his = real(U_his);
LF_his = real(LF_his);

%% Visualize simulation
instdof = -(indp(3)*3-2);
interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','eggbox5x5fold',0.0001,LF_his,instdof,[-inf inf -inf inf])
% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','eggbox5x5fold',0.0001)

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
