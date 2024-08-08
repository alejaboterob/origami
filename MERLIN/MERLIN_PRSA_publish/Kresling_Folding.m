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
%% =========== KRESLING MULTI-STABLE FOLDING =========================== %%
clear all; close all; clc;
%% Define geomtry and material
N = 8; h = 1; lyr = 2; phi = 2*pi/8;
% N-gon prism; height of each layer; number of inter-layer planes; twisting
% angle of the prism in each layer
MaxIcr = 180; blam = 0.032; 
% Maximum increment number & initial load factor
Kf = 1e-3; Kb = Kf; E0 = 5e3; Abar = 0.1;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 45; limrht = 315;
% Left and right limits for the linear range of rotational stiffness
[Node, Panel] = GetDiSym(N,h,lyr,phi);
BarMater = @(Ex)Ogden(Ex, E0); % Define bar material constitutive
RotSpring = @(he,h0,kpi,L0)EnhancedLinear(he,h0,kpi,L0,limlft,limrht);
% Define rotational spring constitutive

%% Set boundary condition
indsupp = find(Node(:,3)<0.01);
nsupp = numel(indsupp);
Supp = [          indsupp(1), 1, 1, 1;
                  indsupp(2), 1, 1, 1;
        indsupp(3:end), zeros(nsupp-2,1)+1, zeros(nsupp-2,1)+1, ones(nsupp-2,1);];
m = size(Node,1);
indp = find(abs(Node(:,3)-max(Node(:,3)))<1e-5); %indp = indp(3);
npp = numel(indp);
Load = [indp, 0*ones(npp,1), 0*ones(npp,1), -1*ones(npp,1);];

%% Perform analysis
[truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarMater,RotSpring,Kf,Kb,Abar);
truss.U0 = zeros(3*size(truss.Node,1),1);
[U_his,LF_his,Data] = PathAnalysis(truss,angles,F,blam,MaxIcr);
U_his = real(U_his);
LF_his = real(LF_his);

%% Visualize simulation
instdof = -indp(1)*3;
interv = 1; endicrm = size(U_his,2);
VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','kreslingfold',0.0001,LF_his,instdof,[-inf inf -inf inf])
% If do not need load-displacement diagram:
% VisualFold(U_his(:,1:interv:endicrm),truss,angles,'none','kreslingfold',0.0001)

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
