% clear all
kpi = 0.0008; h0 = 1*pi;
p = linspace(0.12,2*pi-0.12,1000);
pr = p*0; pk = p*0; pe = p*0;
for i = 1:numel(p)
    [r,k,e] = Spring(p(i),h0,kpi,1);
    pr(i) = r; pk(i) = k; pe(i) = e;
end
% figure(1)
% plot(p,pr)
% ax = gca;
% ax.FontSize = 14;
% ax.XTick = [0 0.5*pi pi 1.5*pi 2*pi];
% ax.XTickLabel = {'0','0.5\pi','\pi','1.5\pi','2\pi'};
% ylim([-inf,inf]); xlim([0,2*pi]); 
% grid on
% hold on
% 
% figure(2)
% plot(p,pk);
% ax = gca;
% ax.FontSize = 14;
% ax.XTick = [0 0.5*pi pi 1.5*pi 2*pi];
% ax.XTickLabel = {'0','0.5\pi','\pi','1.5\pi','2\pi'};
% ylim([-inf,inf]); xlim([0,2*pi]); 
% grid on
% hold on

figure(3)
plot(p,pr)
ax = gca;
ax.FontSize = 14;
ax.XTick = [0 0.5*pi pi 1.5*pi 2*pi];
ax.XTickLabel = {'0','0.5\pi','\pi','1.5\pi','2\pi'};
ylim([-inf,inf]); xlim([0,2*pi]); 
grid on
hold on
