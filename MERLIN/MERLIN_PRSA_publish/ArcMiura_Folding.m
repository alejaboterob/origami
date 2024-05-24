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
%% =========== ARC-MIURA FOLDING ======================================= %%
clear all; close all; clc;
%% Define geomtry and material
load('ArcMiura_45x75')
Node = Node_unfold;
Panel = Panel_unfold;
% Geometry of the Arc_Miura is obtained based on:
% Gattas, J. M., Wu, W., & You, Z. (2013). Miura-base rigid origami: 
% parameterizations of first-level derivative and piecewise geometries. 
% Journal of mechanical design, 135(11), 111011.
MaxIcr = 100; blam = 0.04; 
Kf = 1e-1; Kb = Kf*1e3; E0 = 1e4;  Abar = 1e-1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness
limlft = 45; limrht = 315;
% Left and right limits for the linear range of rotational stiffness
BarMater = @(Ex)Ogden(Ex, E0);  % Define bar material constitutive
RotSpring = @(he,h0,kpi,L0)EnhancedLinear(he,h0,kpi,L0,limlft,limrht);
% Define rotational spring constitutive

%% Set up boundary condition
indsupp = find(Node(:,2)<0.01);
nss = numel(indsupp);
Supp = [          indsupp(1), 1, 1, 1;
                  indsupp(2), 0, 1, 1;
        indsupp(3:end), 0*ones(nss-2,1), ones(nss-2,1), 0*ones(nss-2,1);];
    
m = size(Node,1);
indp = find(abs(Node(:,2)-max(Node(:,2)))<1e-5); 
npp = numel(indp);
Load = [indp, 0*ones(npp,1), -1*ones(npp,1), 0*ones(npp,1);];

%% Perform analysis
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
truss.U0 = zeros(3*size(truss.Node,1),1);
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr);
U_his = real(U_his);
LF_his = real(LF_his);

%% Visualize simulation
instdof = -(indp(1)*3-1);
interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','arcmiurafold',0.0001,LF_his,instdof,[-inf inf -inf inf])
% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','arcmiurafold',0.0001)

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
