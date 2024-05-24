function [Ex,Rbe,Kbe] = BarKe(u,B,L,CM,A)
du = u(1:3)-u(4:6);
Du = [du;-du];
Ex = (B*u/L+0.5*(du'*du)/L^2);
[Sx,Et] = CM(Ex);  
Fx = Sx*A;
if nargout>1, Rbe = Fx*(B'+Du/L); end;
if nargout>2
    Kel = B'*B;
    Kg = Fx/L*[eye(3),-eye(3);-eye(3),eye(3)];
    K1 = ((Du*B)+(Du*B)')/L;
    K2 = (Du*Du')/L^2;
    Kbe = Et*A/L*(Kel+K1+K2)+Kg;
%     Kbe = HA/L*(Kel)+Kg;
end
