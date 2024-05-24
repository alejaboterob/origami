function STAT = PostProcess(Data,truss,angles)
[Sx_bar, ~, Wb] = truss.CM(Data.Exbar);
Rspr_fd = zeros(size(Data.FdAngle)); Efold = Rspr_fd;
Rspr_bd = zeros(size(Data.BdAngle)); Ebend = Rspr_bd;
for i = 1:size(Data.FdAngle,2)
    [Rspr_fdi, ~, Efoldi] = angles.CM(Data.FdAngle(:,i),angles.pf0',angles.kpf',truss.L((size(angles.bend,1)+1):(size(angles.bend,1)+size(angles.fold,1))));
    [Rspr_bdi, ~, Ebendi] = angles.CM(Data.BdAngle(:,i),angles.pb0',angles.kpb',truss.L(1:size(angles.bend,1)));
    Rspr_fd(:,i) = Rspr_fdi; Efold(:,i) = Efoldi;
    Rspr_bd(:,i) = Rspr_bdi; Ebend(:,i) = Ebendi;
end

STAT.bar.Sx = Sx_bar;    
STAT.bar.W = diag(truss.L.*truss.A)*Wb;
STAT.bar.PE = sum(STAT.bar.W,1);

STAT.fold.RM = Rspr_fd;
STAT.fold.E = Efold;
STAT.fold.PE = sum(STAT.fold.E,1);

STAT.bend.RM = Rspr_bd;
STAT.bend.E = Ebend;
STAT.bend.PE = sum(STAT.bend.E,1);

STAT.PE.strain = STAT.bar.PE+STAT.fold.PE+STAT.bend.PE;

