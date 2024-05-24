function [fold, bdry, Trigl] = findfdbd(Panel,bend)
Nn = max(cellfun(@max,Panel)); 
% triangularization
Panelsize = cellfun(@numel,Panel);
Ptri = cell(sum(Panelsize==3),1);
flg = find(Panelsize==3);
for i = 1:sum(Panelsize==3), Ptri{i} = Panel{flg(i)}; end
Triglraw = [bend(:,[1,2,3]);bend(:,[1,2,4]);cell2mat(Ptri)];
Triglraw = sort(Triglraw,2);
Trigl = unique(Triglraw ,'rows');
% formulate connectivity matrix
Comm = sparse(Nn,size(Trigl,1));
for i=1:size(Trigl,1), Comm(Trigl(i,:),i) = true; end;
% search for fold lines
  Ge = Comm'*Comm;
[mf, me] = find(triu(Ge==2)); % triangular meshes that share two common nodes
fold = zeros(length(mf),4);
for i=1:length(mf)
    [link,ia,ib] = intersect(Trigl(mf(i),:),Trigl(me(i),:));
    oftpa = setdiff(1:3,ia);
    oftpb = setdiff(1:3,ib);
    fold(i,:) = [link,Trigl(mf(i),oftpa),Trigl(me(i),oftpb)];
end
fdandbd = sort(fold(:,1:2),2);
onlybd = sort(bend(:,1:2),2);
[~,ibd] = intersect(fdandbd,onlybd,'rows');
fold(ibd,:) = [];

% search for boundaries
Edge = sort([Trigl(:,1) Trigl(:,2); Trigl(:,2) Trigl(:,3); Trigl(:,3) Trigl(:,1)],2);
[u,~,n] = unique(Edge ,'rows');
counts = accumarray(n(:), 1);
bdry = u(counts==1,:);