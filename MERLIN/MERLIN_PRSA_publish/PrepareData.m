function [truss, angles, F] = PrepareData(Node,Panel,Supp,Load,BarCM,RotSpring,kpf,kpb,Abar)
[Bend] = findbend(Panel, Node);
[Fold, Bdry, Trigl] = findfdbd(Panel,Bend);
Bars = [Bend(:,1:2);Fold(:,1:2);Bdry];
[B, L] = dirc3d(Node,Bars);
if size(Supp,1) == 0
    rs = []; 
else
    rs = [reshape([Supp(:,1)*3-2,Supp(:,1)*3-1,Supp(:,1)*3]',[],1),...
          reshape(Supp(:,2:4)',[],1)];
    rs(rs(:,2)==0,:)=[]; rs = rs(:,1);
end

if numel(Abar)==1
    Abar = Abar*ones(size(Bars,1),1);
end

pf0 = zeros(size(Fold,1),1); 
for i = 1:size(Fold,1), pf0(i) = FoldKe(Node,Fold(i,:),kpf,0); end;

pb0 = zeros(size(Bend,1),1); 
for i = 1:size(Bend,1), pb0(i) = FoldKe(Node,Bend(i,:),kpb,0); end;

m = size(Node,1);
F = zeros(3*m,1);
indp = Load(:,1);
F(3*indp-2) = Load(:,2); 
F(3*indp-1) = Load(:,3); 
F(3*indp) = Load(:,4);

truss.CM = BarCM;
truss.Node = Node;
truss.Bars = Bars;
truss.Trigl = Trigl;
truss.B = B; 
truss.L = L;
truss.FixedDofs = unique(rs);
truss.A = Abar;
angles.CM = RotSpring;
angles.fold = Fold;
angles.bend = Bend;
angles.kpf = kpf*ones(1,size(Fold,1));
angles.kpb = kpb*ones(1,size(Bend,1));
angles.pf0 = pf0';
angles.pb0 = pb0'*0+pi;
angles.Panel = Panel;

