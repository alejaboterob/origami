function [Rspr, Kspr, Espr] = EnhancedLinear(he,h0,kpi,L0,limlft,limrht)
limlft = limlft/180*pi; partl = pi/limlft;
limrht = limrht/180*pi; partr = pi/(2*pi-limrht);
% limlft: theta_1: left partition point
% limrht: theta_2: right partition point
if numel(kpi)==1, kpi = kpi*ones(size(he)); end;
Rspr = zeros(size(he)); Kspr = Rspr; 

Lind = he<limlft; Rind = he>limrht; Mind = ~(Lind|Rind);
Rspr(Lind) = kpi(Lind).*real(limlft-h0(Lind))+kpi(Lind).*tan(partl/2*(he(Lind)-limlft))/(partl/2);
Kspr(Lind) = kpi(Lind).*sec(partl/2*(he(Lind)-limlft)).^2;
Rspr(Rind) = kpi(Rind).*real(limrht-h0(Rind))+kpi(Rind).*tan(partr/2*(he(Rind)-limrht))/(partr/2);
Kspr(Rind) = kpi(Rind).*sec(partr/2*(he(Rind)-limrht)).^2;
Rspr(Mind) = kpi(Mind).*real(he(Mind)-h0(Mind));
Kspr(Mind) = kpi(Mind);
Rspr = L0.*Rspr; Kspr = L0.*Kspr;

if nargout>2
    Espr = zeros(size(he));
    Espr(Lind) = 0.5*kpi(Lind).*real(h0(Lind)-limlft).^2+kpi(Lind).*real(h0(Lind)-limlft).*(limlft-he(Lind))-4*kpi(Lind)/partl^2.*log(abs(cos(partl/2*(limlft-he(Lind)))));
    Espr(Rind) = 0.5*kpi(Rind).*real(limrht-h0(Rind)).^2+kpi(Rind).*real(limrht-h0(Rind)).*(he(Rind)-limrht)-4*kpi(Rind)/partr^2.*log(abs(cos(partr/2*(he(Rind)-limrht))));
    Espr(Mind) = 0.5*kpi(Mind).*real(he(Mind)-h0(Mind)).^2;
    Espr = L0.*Espr;
end
