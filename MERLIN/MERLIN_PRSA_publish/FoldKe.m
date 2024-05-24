function [he,Rhe,Khe] = FoldKe(Cood, List, kpi, h0, L0, CM)
rkj = [Cood(List(2),:)-Cood(List(1),:)]'; 
rij = [Cood(List(3),:)-Cood(List(1),:)]'; 
rkl = [Cood(List(2),:)-Cood(List(4),:)]'; 
rmj = icross(rij,rkj); rnk = icross(rkj,rkl);
sgn = ((abs(rnk'*rij)>1e-8)*sign(rnk'*rij)+(abs(rnk'*rij)<=1e-8)*1);
he = real(acos(rmj'*rnk/(norm(rmj)*norm(rnk)))); 
he = real(sgn*he);
if he<0 
    he = 2*pi+he; 
end;

if nargin > 4
[Rspr, Kspr] = CM(he,h0,kpi,L0);

if nargout>1
    di = norm(rkj)/(rmj'*rmj)*rmj;
    dl = -norm(rkj)/(rnk'*rnk)*rnk;
    dj = (rij'*rkj/(rkj'*rkj)-1)*di-rkl'*rkj/(rkj'*rkj)*dl;
    dk = -rij'*rkj/(rkj'*rkj)*di+(rkl'*rkj/(rkj'*rkj)-1)*dl;

    Jhe = [dj;dk;di;dl];
    Rhe = Rspr*Jhe;
end

if nargout > 2
    dii = -norm(rkj)/(rmj'*rmj)^2*((rmj*cross(rkj,rmj)')+(rmj*cross(rkj,rmj)')');
    
    dtempij = -norm(rkj)/(rmj'*rmj)^2*(rmj*(cross(rij-rkj,rmj))'+(cross(rij-rkj,rmj))*rmj');
    dij = -rmj*rkj'/(rmj'*rmj*norm(rkj))+dtempij;
    
    dtempik = norm(rkj)/(rmj'*rmj)^2*(rmj*(cross(rij,rmj))'+(cross(rij,rmj))*rmj');
    dik = rmj*rkj'/(rmj'*rmj*norm(rkj))+dtempik;
    
    dil = zeros(3);
    
    dll = norm(rkj)/(rnk'*rnk)^2*(rnk*cross(rkj,rnk)'+(rnk*cross(rkj,rnk)')');
    
    dtemplk = norm(rkj)/(rnk'*rnk)^2*(rnk*(cross(rkl-rkj,rnk))'+(cross(rkl-rkj,rnk))*rnk');
    dlk = -rnk*rkj'/(rnk'*rnk*norm(rkj))+dtemplk;
    
    dtemplj = norm(rkj)/(rnk'*rnk)^2*(rnk*(cross(rnk,rkl))'+(rnk*(cross(rnk,rkl))')');
    dlj = rnk*rkj'/(rnk'*rnk*norm(rkj))+dtemplj;
    
    dT1jj = 1/(rkj'*rkj)*((-1+2*rij'*rkj/(rkj'*rkj))*rkj-rij);
    dT2jj = 1/(rkj'*rkj)*(2*rkl'*rkj/(rkj'*rkj)*rkj-rkl);
    djj = di*dT1jj'+(rij'*rkj/(rkj'*rkj)-1)*dij-(dl*dT2jj'+rkl'*rkj/(rkj'*rkj)*dlj);
    
    dT1jk = 1/(rkj'*rkj)*(-2*rij'*rkj/(rkj'*rkj)*rkj+rij);
    dT2jk = 1/(rkj'*rkj)*((1-2*rkl'*rkj/(rkj'*rkj))*rkj+rkl);
    djk = di*dT1jk'+(rij'*rkj/(rkj'*rkj)-1)*dik-(dl*dT2jk'+rkl'*rkj/(rkj'*rkj)*dlk);
    
    dT1kk = dT2jk;
    dT2kk = dT1jk;
    dkk = dl*dT1kk'+(rkl'*rkj/(rkj'*rkj)-1)*dlk-(di*dT2kk'+rij'*rkj/(rkj'*rkj)*dik);
    
    Hp = [ djj , djk , dij', dlj';
           djk', dkk , dik', dlk';
           dij , dik , dii , dil ;
           dlj , dlk , dil', dll];
          
                                  
    Khe = (Kspr*(Jhe*Jhe')+Rspr*Hp);
end
end