function [Sx, Et, Wb] = Ogden(Ex, C0)
% Ogden hyperelastic constitutive model for bar elements
alfa = [5,1]; % Specify parameteres
pstr = real(sqrt(2*Ex+1));
C0 = (pstr<1)*1*C0+(~(pstr<1))*1*C0;
Et = C0/(alfa(1)-alfa(2)).*((alfa(1)-2)*pstr.^(alfa(1)-4)-(alfa(2)-2)*pstr.^(alfa(2)-4));
Sx = C0/(alfa(1)-alfa(2)).*(pstr.^(alfa(1)-2)-pstr.^(alfa(2)-2));
if nargout>2
    Wb = C0/(alfa(1)-alfa(2)).*((pstr.^alfa(1)-1)/alfa(1)-(pstr.^alfa(2)-1)/alfa(2));
end
    