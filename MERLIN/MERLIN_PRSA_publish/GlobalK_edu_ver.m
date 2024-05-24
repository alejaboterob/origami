function [IF,K] = GlobalK_edu_ver(Ui,Node,truss,angles)
Nn = size(Node,1);
IFb = zeros(3*Nn,1); IFp = IFb;
indi = zeros(36*size(truss.Bars,1),1); indj = indi; kentry = indi;
Nodenw(:,1) = Node(:,1)+Ui(1:3:end);
Nodenw(:,2) = Node(:,2)+Ui(2:3:end);
Nodenw(:,3) = Node(:,3)+Ui(3:3:end);

for bel = 1:size(truss.Bars,1) 
    eDof = [(-2:0)+(truss.Bars(bel,1)*3),(-2:0)+(truss.Bars(bel,2)*3)]';
    [~,Rbe,Kbe] = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel));
    IFb(eDof) = IFb(eDof)+Rbe;
    I=repmat(eDof,1,6); J=I';
    indi(36*(bel-1)+1:36*bel) = I(:);
    indj(36*(bel-1)+1:36*bel) = J(:); 
    kentry(36*(bel-1)+1:36*bel) = Kbe(:);
end
Kb = sparse(indi,indj,kentry,3*Nn,3*Nn);

indi = zeros(144*size(angles.bend,1),1); indj = indi; kentry = indi;
Lbend = truss.L(1:size(angles.bend,1));
for del = 1:size(angles.bend,1)
    eDof = reshape([3*angles.bend(del,:)-2;...
                    3*angles.bend(del,:)-1;...
                    3*angles.bend(del,:)],12,1);
    bend = angles.bend(del,:);
    [~,Rpe,Kpe] = FoldKe(Nodenw,bend,angles.kpb,angles.pb0(del),Lbend(del),angles.CM);
    IFp(eDof) = IFp(eDof)+Rpe;
    I=repmat(eDof,1,12); J=I';
    indi(144*(del-1)+1:144*del) = I(:);
    indj(144*(del-1)+1:144*del) = J(:); 
    kentry(144*(del-1)+1:144*del) = Kpe(:);
end;
Kbd = sparse(indi,indj,kentry,3*Nn,3*Nn);
if isempty(Kbd), Kbd = zeros(3*Nn); end

indi = zeros(144*size(angles.fold,1),1); indj = indi; kentry = indi;
Lfold = truss.L(size(angles.bend,1)+1:size(angles.bend,1)+size(angles.fold,1));
for fel = 1:size(angles.fold,1)
    eDof = reshape([3*angles.fold(fel,:)-2;...
                    3*angles.fold(fel,:)-1;...
                    3*angles.fold(fel,:)],12,1);
    fold = angles.fold(fel,:);
    [~,Rpe,Kpe] = FoldKe(Nodenw,fold,angles.kpf,angles.pf0(fel),Lfold(fel),angles.CM);
    IFp(eDof) = IFp(eDof)+Rpe;
    I=repmat(eDof,1,12); J=I';
    indi(144*(fel-1)+1:144*fel) = I(:);
    indj(144*(fel-1)+1:144*fel) = J(:); 
    kentry(144*(fel-1)+1:144*fel) = Kpe(:);
end;
Kfd = sparse(indi,indj,kentry,3*Nn,3*Nn);

IF = IFb+IFp;
K = Kb+Kbd+Kfd;
K = (K+K')/2;