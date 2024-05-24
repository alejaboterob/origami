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
%% =========== MIURA POP-THROUGH ======================================= %%
clear all; close all; clc;
%% Define geomtry and material
sec_hor=5;  sec_vert=5; 
theta = 60; a = 2; b = 2; fdang = 30; 
% Geometry of the Miura: a, b are the edge lengths of each parallelogram
% panel; fdang controls the folding angle; theta is the panel angle
MaxIcr = 160; blam = 0.06; 
% Maximum increment number & initial load factor
Kf = 1e-1; Kb = Kf*10; E0 = 1e4; Abar = 1e-1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 45; limrht = 315;
% Left and right limits for the linear range of rotational stiffness
[Node,Panel,~]=ConfigMiura(sec_hor,sec_vert,theta,a,b,fdang); 
BarMater = @(Ex)Ogden(Ex, E0); % Define bar material constitutive
RotSpring = @(he,h0,kpi,L0)EnhancedLinear(he,h0,kpi,L0,limlft,limrht);

%% Set up boundary conditions
m = size(Node,1);
leftxz = [3:2:(sec_vert*2+1)]';
leftx = [2:2:(sec_vert*2+1)]';
rightz = [1;leftxz]+(sec_vert*2+1)*(sec_hor*2);
Supp = [          1, 1, 1, 1; 
        leftxz, 0*leftxz+1, 0*leftxz, 0*leftxz+1;
        rightz(1), 0, 1, 1; 
        rightz(2:end), 0*rightz(2:end), 0*rightz(2:end), 0*rightz(2:end)+1];
indp = 61;
ff = -1*ones(length(indp),1); 
Load = [indp,zeros(length(indp),1),zeros(length(indp),1),ff];
indp = Load(:,1);

%% Perform analysis
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
truss.U0 = zeros(3*size(truss.Node,1),1);
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr);
U_his = real(U_his);
LF_his = real(LF_his);

%% Visualize simulation
instdof = -indp(1)*3;
interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','miura5x5ptd',0.0001,LF_his,instdof,[-inf inf -inf inf])
% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','miura5x5ptd',0.0001)

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
