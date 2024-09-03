clc
clear
%��������
N=1000;
lamda1=0.8;
lamda2=0.2;
delta1=0.5;
delta2=0.5;
mu=0.5;
gama1=2;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
gama2=0.5;
% ma=5;%AD����ÿ�����ӵı���
mb=5;
ter=20;% mu��ȡֵ����
stp=100;%ʱ�䳤��
mu=0.5;
yita=1;%��������

Exponent_B=2.1;%AD�����Ծָ��
Exponent_A=2.1;%AD�����Ծָ��
betaR_termi=40;
ma_termi=20;
 %��������������󣬱������仯
 load AActivity ;
 load BActivity ;
  rhoG=zeros(betaR_termi,ma_termi);
parfor la=1:betaR_termi  
beta_R=la/(2*betaR_termi);
  for xun=1:ma_termi      
 ma=xun/2;
%ÿ��Exponent_B�����������������������Ҫ��η���
%������Ծ�Ⱦ���ÿһ��gamma��������µĻ�Ծ�Ⱦ��󣬲���ÿһ�η��涼�����������������Ӧ����ͬһ��gamma��ͬ����������
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
            %�����ڽӾ���
           
%����������
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
        
          rhoG(la,xun)=sum(PG)/N;
end
 end

x=(1:40)/2;
y=(1:40)/(2*40);
image(x,y,rhoG,'CDataMapping','scaled');
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$m_a$');
ylabel('$\beta$');
 saveas(gcf, 'MMCA.fig'); %���浱ǰ���ڵ�ͼ��
 save('rhoG.mat','rhoG')