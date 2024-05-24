%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 MERLIN2                               %%
%                     Written by: Ke Liu (ke.liu@gatech.edu)              %
% Ref: K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid    %
%      origami - An efficient computational approach.' PRSA.              %
%      K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to  % 
%      capture highly nonlinear behavior of non-rigid origami.'           %
%      Proceedings of IASS Annual Symposium 2016.                         %
%      E. T. Filipov, K. Liu, T. Tachi, M. Schenk, G. H. Paulino (2017).  %
%      'Bar and hinge models for scalable analysis of origami.'  IJSS     %
%      K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear          %
%      structural analysis of origami assemblages using the MERLIN2       %
%      software.' Origami^7.                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =========== MIURA BENDING =========================================== %%
clear all; % close all; clc;
%% Define geomtry
Node = [ 0,  0,  0;
         0, 10,  0;
        12, 15,  0;
        30, 10,  0;
        30,  0,  0;
        18, -5,  0]*2;
Panel = {[1:6]};

%% Set up boundary conditions
m = size(Node,1);
Supp = [ 1, 1, 1, 1;
         2, 1, 0, 1;
         4, 0, 0, 1];
Load = [ 5, 0, 0,  -1];

%% Adopt generalized N5B8 model
AnalyInputOpt = struct(...
    'ModelType','N5B8',...
    'MaterCalib','auto',...  % 'manual'
    'ModElastic', 5e3,...
    'Poisson', 0.35,...
    'Thickness', 0.127,... 
    'LScaleFactor', 3,...
    'LoadType','Force',...
    'InitialLoadFactor', 0.00001,...
    'MaxIcr', 100,...
    'StopCriterion',@(Node,U,icrm)(abs(U(5*3))>12));

%% Adopt generalized N4B5 model
% AnalyInputOpt = struct(...
%     'ModelType','N4B5',...
%     'MaterCalib','manual',...  % 'manual'
%     'BarCM', @(Ex)Ogden(Ex, 5e3),...
%     'Abar', 2.4652,...
%     'Kb',0.9726*([38.4187;38.4187;41.7612]/0.127).^(1/3)./[38.4187;38.4187;41.7612],...
%     'Kf',1,...
%     'RotSprBend', @SuperLinearBend,...
%     'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,15,345),...
%     'LoadType','Force',...
%     'InitialLoadFactor', 0.00001,...
%     'MaxIcr', 100,...
%     'StopCriterion',@(Node,U,icrm)(abs(U(5*3))>12));

%% Perform analysis
% Assemble input data
[truss, angles, AnalyInputOpt] = PrepareData(Node,Panel,Supp,Load,AnalyInputOpt);
% Specify initial deformation state
truss.U0 = zeros(3*size(truss.Node,1),1);
% Perform path-following analysis using MGDCM
[Uhis,Fhis] = PathAnalysis(truss,angles,AnalyInputOpt);
% Clean output data
Uhis = real(Uhis);
Fhis = real(Fhis);
STAT = PostProcess(Uhis,truss,angles); 

%% Visualize simulation
instdof = [5,-3];
interv = 1; endicrm = size(Uhis,2);
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,Fhis(1:interv:endicrm,:),instdof)
% If do not need load-displacement diagram:
% VisualFold(Uhis(:,1:interv:endicrm),truss,angles,[],[])
% To visualize the strains in bars:
% VisualFold(Uhis(:,1:10:endicrm),truss,angles,[],[],'IntensityMap','Edge','IntensityData',STAT.bar.Sx(:,1:10:endicrm),'ShowInitial','off')

%% Plot diagrams
% Force vs Displacement
figure()
dsp = sign(instdof(2))*Uhis((instdof(1)*3-(3-abs(instdof(2)))),:);
plot(dsp,Fhis)
hold on
axis tight

% Strain energy vs Displacement
figure()
plot(dsp,STAT.PE,'r-','linewidth',2);
axis tight
xlabel('Displacement','fontsize',14);
ylabel('Stored Energy','fontsize',14);

%% Plot final configuration
Ux = Uhis(:,end);
Nodew = truss.Node;
Nodew(:,1) = truss.Node(:,1)+Ux(1:3:end);
Nodew(:,2) = truss.Node(:,2)+Ux(2:3:end);
Nodew(:,3) = truss.Node(:,3)+Ux(3:3:end); 
figure()
PlotOri(truss.Node,angles.Panel,truss.Trigl,'FoldEdgeStyle','-','EdgeShade',0.3,'PanelColor','none');
PlotOri(Nodew,angles.Panel,truss.Trigl,'PanelColor','g');
axis equal; axis off;
camproj('perspective')
light
view(-2,16)
rotate3d on
