clear all
E0 = 1; 
%% Neo-Hookean
% alfa = [2,0];
%% St. Venant
% alfa = [4,2];
%% Others
% alfa = [5,1];

pstr = (0.1:0.01:1.9)';
Ex = 0.5*(pstr.^2-1);
[Sx, Et, Wb] = Ogden(Ex,E0);
% Sx = @(pstr,E0)E0/(alfa(1)-alfa(2))*(pstr.^(alfa(1)-2)-pstr.^(alfa(2)-2));
% Wb = @(pstr,E0)E0/(alfa(1)-alfa(2))*((pstr.^alfa(1)-1)/alfa(1)-(pstr.^alfa(2)-1)/alfa(2));
% Wb = @(pstr,E0)E0/(alfa(1)-alfa(2))*(pstr.^alfa(1)-1)/alfa(1);

figure(21)
hold on
plot(pstr,Sx,'linewidth',2);
xlabel('\lambda','fontsize',14)
ylabel('S_x_x','fontsize',14)
axis([0.7 1.3 -inf inf]);
hold on