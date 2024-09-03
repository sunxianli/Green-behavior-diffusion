clc
clear
%deta1改变
N=1000;
lamda1=0.6;
lamda2=0.6;
% delta1=0.3;
delta2=0.3;
med=[0.1,0.3,0.5,0.7,0.9];
mu=0.5;
gama1=1.5;
gama2=0.5;
ma=5;
mb=5;
ter=40;%beta的取值个数
stp=150;
rhoG=zeros(length(med),ter);
load('BActivity.mat')
load('AActivity.mat')
parfor m1=1:length(med);
   delta1=med(m1);
            for xun = 1:ter
    beta_R=xun/80;

            PA2S=0.01*ones(1,N);
            PUS=0.98*ones(1,N);
            PA1G=0.01*ones(1,N);
            PA1S=0*ones(1,N);

            PUS_UPDATE=zeros(1,N);
            PA2S_UPDATE=zeros(1,N);
            PA1G_UPDATE=zeros(1,N);
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);
            rA2=zeros(1,N);
            qA1=zeros(1,N);
            qA2=zeros(1,N);
            qU=zeros(1,N);

            RA1=zeros(N,N);
            RA2=zeros(N,N);
            QA1=zeros(N,N);
            QA2=zeros(N,N);
            QU=zeros(N,N);

            for t = 1:stp
                for i =1:N
                    for j =1:N
                        RA1(j,i)=1-(AActivity(i)+AActivity(j))*(PA1G(1,j)+PA1S(1,j))*lamda1*ma/N;
                        RA2(j,i)=1-(AActivity(i)+AActivity(j))*PA2S(1,j)*lamda2*ma/N;
                        QA1(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*gama1*mb/N;
                        QA2(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*gama2*mb/N;
                        QU(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*mb/N;
                    end
                    tempprodRA1=cumprod(RA1(:,i));
                    rA1(1,i)=tempprodRA1(N);
                    tempprodRA2=cumprod(RA2(:,i));
                    rA2(1,i)=tempprodRA2(N);
                    tempprodQA1=cumprod(QA1(:,i));
                    qA1(1,i)=tempprodQA1(N);
                    tempprodQA2=cumprod(QA2(:,i));
                    qA2(1,i)=tempprodQA2(N); 
                    tempprodQU=cumprod(QU(:,i));
                    qU(1,i)=tempprodQU(N);
 PUS_UPDATE(1,i)=PA1S(1,i)*delta1*qU(1,i)+PA2S(1,i)*delta2*qU(1,i)+PA1G(1,i)*delta1*mu+PUS(1,i)* rA1(1,i)*rA2(1,i)*qU(1,i);
 PA1G_UPDATE(1,i)=PA1S(1,i)*(delta1*(1-qU(1,i))+(1-delta1)*(1-qA1(1,i)))+PA2S(1,i)*(delta2*(1-qU(1,i))+(1-delta2)*(1-qA2(1,i)))...
                 +PUS(1,i)*(rA1(1,i)*rA2(1,i)*(1-qU(1,i))+(1-rA1(1,i)-(1-rA1(1,i))*(1-rA2(1,i)))*(1-qA1(1,i))+(1-rA2(1,i))*(1-qA2(1,i)))...
                 +PA1G(1,i)*(1-mu);
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1)*qA1(1,i)+PA1G(1,i)*(1-delta1)*mu+PUS(1,i)*(1-rA1(1,i)-(1-rA1(1,i))*(1-rA2(1,i)))*qA1(1,i);  
 PA2S_UPDATE(1,i)=PA2S(1,i)*(1-delta2)*qA2(1,i)+PUS(1,i)*(1-rA2(1,i))*qA2(1,i); 

                end
                PUS=PUS_UPDATE;
                PA1G=PA1G_UPDATE;
                PA1S= PA1S_UPDATE;
                PA2S=PA2S_UPDATE;
            end
 rhoG(m1,xun)=sum(PA1G)/N;

            end  
end
xzhou=(1:ter)/80;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
set(gca, 'TickLabelInterpreter', 'latex');
plot(xzhou,rhoG(1,:)','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,rhoG(2,:)','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
plot(xzhou,rhoG(3,:),'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
plot(xzhou,rhoG(4,:)','-d','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,rhoG(5,:),'-+','color',[238/256 228/256 132/256],'MarkerFaceColor',[238/256 228/256 132/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',16.5);ylabel('$\rho^{G}$','FontSize',16.5);   % é§æ­ç£æç£Ðé²ï¿?
h=legend('$\mu=0.1$','$\mu=0.3$','$\mu=0.5$','$\mu=0.7$','$\mu=0.9$');
set(h,'Interpreter','latex','FontSize',16.5)%,
save('rhoG.mat','rhoG')
saveas(gcf, 'save.fig')

