ter=20;
load rhoG_ave.mat
xzhou=1:ter;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
plot(xzhou,rhoG_num(1,:));
plot(xzhou,rhoG_num(2,:));
plot(xzhou,rhoG_num(3,:));
plot(xzhou,rhoG_num(4,:));
plot(xzhou,rhoG_num(5,:));
plot(xzhou,rhoG_num(6,:));
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$m_p$','FontSize',15);ylabel('$\rho^G$','FontSize',15);   % 鍧愭爣杞磋В閲�?
h=legend('$\lambda_1=0.9,\lambda_2=0.1$','$\lambda_1=0.8,\lambda_2=0.2$','$\lambda_1=0.7,\lambda_2=0.3$','$\lambda_1=0.1,\lambda_2=0.9$','$\lambda_1=0.2,\lambda_2=0.8$','$\lambda_1=0.3,\lambda_2=0.7$');
set(h,'Interpreter','latex','FontSize',15)%,
save('rhoG_ave.mat','rhoG_num')
saveas(gcf, 'rho.fig')