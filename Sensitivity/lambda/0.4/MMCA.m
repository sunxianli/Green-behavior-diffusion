clc
clear
%参数设置
N=1000;
delta1=0.3;
delta2=0.3;
mu=0.5;
gama1=1.2;%衰减系数，当其为0时，知道信息的个体是完全免疫的。在这里若等于1，
gama2=0.2;
ma=5;%AD网络每次连接的边数
mb=5;
stp=150;%时间长度
%定义矩阵
lamda1_termi=40;
lamda2_termi=40;
rhoG=zeros(lamda1_termi,lamda2_termi);;%节点G的密度
beta_R=0.4;
load AActivity 
load BActivity 
  parfor l = 1:lamda2_termi
     lamda2= l/40;
     for xun= 1:lamda1_termi
     lamda1= xun/40;
     
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
%计算马氏链
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
%             PA1=PA1S+PA1G;
%             PA2=PA2S;

           rhoG(l,xun)=sum(PG)/N;

        end
  end
save('rhoG.mat','rhoG');
x=(1:40)/40;
y=(1:40)/40;
image(x,y,rhoG,'CDataMapping','scaled');
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\lambda_1$');
ylabel('$\lambda_2$');
saveas(gcf, 'MMCA.fig'); %保存当前窗口的图像
