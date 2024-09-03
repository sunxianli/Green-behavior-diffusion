load rhoG.mat
ter=40;
xzhou=(1:ter)/80;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
set(gca, 'TickLabelInterpreter', 'latex');
plot(xzhou,rhoG(1,:),'-+','color',[0.64 0.08 0.18],'MarkerFaceColor',[0.64 0.08 0.18]);
plot(xzhou,rhoG(2,:)','-d','color',[180/256 139/256 171/256],'MarkerFaceColor',[180/256 139/256 171/256]);
plot(xzhou,rhoG(3,:)','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,rhoG(4,:)','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
plot(xzhou,rhoG(5,:),'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);

set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',16.5);ylabel('$\rho^{G}$','FontSize',16.5);   % é§æ­ç£æç£Ðé²ï¿?
h=legend('$\mu=0.1$','$\mu=0.3$','$\mu=0.5$','$\mu=0.7$','$\mu=0.9$');
set(h,'Interpreter','latex','FontSize',16.5)%,
