clc
clear
%鍙傛暟璁剧疆
N=1000;
lamda1=0.6;
lamda2=0.6;
delta1=0.3;
delta2=0.3;
mu=0.5;
gama1=1.5;%琛板噺绯绘暟锛屽綋鍏朵负0鏃讹紝鐭ラ亾淇℃伅鐨勪釜浣撴槸�?�屽叏鍏嶇柅鐨勩�傚湪杩欓噷鑻ョ瓑浜�?1锛�
gama2=0.5;
ma=5;%AD缃戠粶姣忔杩炴帴鐨勮竟鏁�
mb=5;
ter=40;% beta鐨勫彇鍊间釜鏁�
MMCA_rep=1;%浠跨湡娆℃暟
stp=150;%鏃堕棿闀垮害
%瀹氫箟鐭╅樀
rhoG=zeros(1,ter);%鑺傜偣G鐨勫瘑搴�?
rhoA1=zeros(1,ter);
rhoA2=zeros(1,ter);
load('BActivity.mat')
load('AActivity.mat')
            parfor xun = 1:ter
    beta_R=xun/80;
        for rep = 1:MMCA_rep
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
            %鐢熸垚閭绘帴鐭╅�?
           
%璁＄畻椹皬閾�
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

            PG=PA1G;
            PA1=PA1S+PA1G;
            PA2=PA2S;

          rhoG(xun)=rhoG(xun)+sum(PG)/N/MMCA_rep;
          rhoA1(xun)=rhoA1(xun)+sum(PA1)/N/MMCA_rep;
          rhoA2(xun)=rhoA2(xun)+sum(PA2)/N/MMCA_rep;
        end
end

xzhou=(1:ter)/80;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
plot(xzhou,rhoG','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,rhoA1','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
plot(xzhou,rhoA2,'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',15);ylabel('proportion','FontSize',15);   % 鍧愭爣杞磋В閲�?
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,
save('rhoG.mat','rhoG')
save('rhoA1.mat','rhoA1')
save('rhoA2.mat','rhoA2')
saveas(gcf, 'save.fig')

