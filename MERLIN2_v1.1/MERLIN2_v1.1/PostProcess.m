function STAT = PostProcess(U_his,truss,angles)
%% Get Data
Exbar = zeros(size(truss.Bars,1),size(U_his,2)); 
FdAngle = zeros(size(angles.fold,1),size(U_his,2)); 
BdAngle = zeros(size(angles.bend,1),size(U_his,2));
for icrm=1:size(U_his,2)
    Ui = U_his(:,icrm);
    Nodenw = truss.Node;
    Nodenw(:,1) = truss.Node(:,1)+Ui(1:3:end);
    Nodenw(:,2) = truss.Node(:,2)+Ui(2:3:end);
    Nodenw(:,3) = truss.Node(:,3)+Ui(3:3:end);
    
    eDofb = kron(truss.Bars,3*ones(1,3))+repmat([-2,-1,0],size(truss.Bars,1),2);
    du = [Ui(eDofb(:,1:3))-Ui(eDofb(:,4:6))];
    Exbar(:,icrm) = truss.B*Ui./truss.L+0.5*sum(du.^2,2)./(truss.L.^2);

    for del = 1:size(angles.bend,1)
        bend = angles.bend(del,:);
        BdAngle(del,icrm) = FoldKe(Nodenw,bend);
    end

    for fel = 1:size(angles.fold,1)
        fold = angles.fold(fel,:);
        FdAngle(fel,icrm) = FoldKe(Nodenw,fold);
    end
end

%% Interpret Data
[Sx_bar, ~, Wb] = truss.CM(Exbar);
Rspr_fd = zeros(size(FdAngle)); Efold = Rspr_fd;
Rspr_bd = zeros(size(BdAngle)); Ebend = Rspr_bd;
for i = 1:size(U_his,2)
    [Rspr_fdi, ~, Efoldi] = angles.CMfold(FdAngle(:,i),angles.pf0,angles.Kf,truss.L((size(angles.bend,1)+1):(size(angles.bend,1)+size(angles.fold,1))));
    [Rspr_bdi, ~, Ebendi] = angles.CMbend(BdAngle(:,i),angles.pb0,angles.Kb,truss.L(1:size(angles.bend,1)));
    Rspr_fd(:,i) = Rspr_fdi; Efold(:,i) = Efoldi;
    Rspr_bd(:,i) = Rspr_bdi; Ebend(:,i) = Ebendi;
end

STAT.bar.Ex = Exbar; 
STAT.bar.Sx = Sx_bar;    
STAT.bar.USi = diag(truss.L.*truss.A)*Wb;
STAT.bar.US = sum(STAT.bar.USi,1);

STAT.fold.Angle = FdAngle;
STAT.fold.RM = Rspr_fd;
STAT.fold.UFi = Efold;
STAT.fold.UF = sum(STAT.fold.UFi,1);

STAT.bend.Angle = BdAngle;
STAT.bend.RM = Rspr_bd;
STAT.bend.UBi = Ebend;
STAT.bend.UB = sum(STAT.bend.UBi,1);

STAT.PE = STAT.bar.US+STAT.fold.UF+STAT.bend.UB;

