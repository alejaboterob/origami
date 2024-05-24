function [Exbar,FdAngle,BdAngle,LFd,LBd] = GetData(Ui,Node,truss,angles)
Exbar = zeros(size(truss.Bars,1),1); 
FdAngle = zeros(size(angles.fold,1),1); LFd = FdAngle;
BdAngle = zeros(size(angles.bend,1),1); LBd = BdAngle;
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+Ui(1:3:end);
Nodenw(:,2) = Node(:,2)+Ui(2:3:end);
Nodenw(:,3) = Node(:,3)+Ui(3:3:end);
for bel = 1:size(truss.Bars,1) 
    eDof = [(-2:0)+(truss.Bars(bel,1)*3),(-2:0)+(truss.Bars(bel,2)*3)]';
    Exbar(bel) = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel));
end

for del = 1:size(angles.bend,1)
    bend = angles.bend(del,:);
    BdAngle(del) = FoldKe(Nodenw,bend,angles.kpb,angles.pb0(del));
    LBd(del) = norm(Nodenw(bend(2),:)-Nodenw(bend(1),:));
end;

for fel = 1:size(angles.fold,1)
    fold = angles.fold(fel,:);
    FdAngle(fel) = FoldKe(Nodenw,fold,angles.kpf,angles.pf0(fel));
    LFd(fel) = norm(Nodenw(fold(2),:)-Nodenw(fold(1),:));
end;