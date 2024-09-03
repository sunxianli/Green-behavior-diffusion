load rhoG.mat
load rhoA1.mat
load rhoA2.mat
ter=40;
xzhou=(1:ter)/80;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
plot(xzhou,rhoG','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,rhoA1','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
plot(xzhou,rhoA2,'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',15);ylabel('proportion','FontSize',15);   % é§æ­ç£æç£Ðé²ï¿?
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,
saveas(gcf, 'save1505.fig')
