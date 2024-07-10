%ʱ���A1A2��һ�η���
clc
clear
%��������
N=1000;
lamda1=0.8;
lamda2=0.2;
delta1=0.5;
delta2=0.5;
gama1=2;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
gama2=0.5;

ma=5;%AD����ÿ�����ӵı���
mb=5;
ter=20;% mu��ȡֵ����
stp=100;%ʱ�䳤��
mu=0.5;
yita=1;%��������
minActivity=0.1;%�½�
% Exponent_A=2.1;%AD�����Ծָ��
 Exponent_B=2.1;
%�������
MMCA_rep=2;
beta=zeros(1,ter);%��ͬgamma�µ���ֵ���
beta_MMCA=zeros(MMCA_rep,ter);
beta_ave=zeros(1,ter);
% beta_num=zeros(num,ter);


for rep = 1:MMCA_rep
%Ԥ�������������
  temAct=rand(1,N);
  temAct2=rand(1,N);
for xun = 1:ter
 Exponent_A=(xun+13)/7;
%ÿ��Exponent_B�����������������������Ҫ��η���
%������Ծ�Ⱦ���ÿһ��gamma��������µĻ�Ծ�Ⱦ��󣬲���ÿһ�η��涼�����������������Ӧ����ͬһ��gamma��ͬ����������
AActivity=1:N;
BActivity=1:N;
   for i = 1:N
      AActivity(i)=yita*power(((power(1,-Exponent_A+1)-power(minActivity,-Exponent_A+1))*temAct(i)+power(minActivity,-Exponent_A+1)),1/(-Exponent_A+1));
   end
   for i = 1:N 
 BActivity(i)=yita*power(((power(1,-Exponent_B+1)-power(minActivity,-Exponent_B+1))*temAct2(i)+power(minActivity,-Exponent_B+1)),1/(-Exponent_B+1));
end
   %���ó�ʼֵ����ʼ����
           
            PA2S=0.01*ones(1,N);
            PA1S=0.01*ones(1,N);

            PA2S_UPDATE=zeros(1,N);
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);
            rA2=zeros(1,N);

            RA1=zeros(N,N);
            RA2=zeros(N,N);

            %�����ڽӾ���
         
%����������
            for t = 1:stp
                for i =1:N
                    for j =1:N 
                        RA1(j,i)=1-(AActivity(i)+AActivity(j))*PA1S(1,j)*lamda1*ma/N;
                        RA2(j,i)=1-(AActivity(i)+AActivity(j))*PA2S(1,j)*lamda2*ma/N;
                    end
                    tempprodRA1=cumprod(RA1(:,i));
                    rA1(1,i)=tempprodRA1(N);
                    tempprodRA2=cumprod(RA2(:,i));
                    rA2(1,i)=tempprodRA2(N);
                    
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1)+(1-PA1S(1,i)-PA2S(1,i))*(rA2(1,i)-rA1(1,i)*rA2(1,i));  
 PA2S_UPDATE(1,i)=PA2S(1,i)*(1-delta2)+(1-PA1S(1,i)-PA2S(1,i))*(1-rA2(1,i)); 
                end
                PA1S= PA1S_UPDATE;
                PA2S=PA2S_UPDATE;
            end       
%����ֵ��ͼ��ʱ����ϵͳ������̬ʱ��G��ֵΪ0.����A1=A1S
  PA1=PA1S;
  PA2=PA2S;
    
 rhoA1=sum(PA1)/N;
 rhoA2=sum(PA2)/N;%�����ĳ�����µľ�ֵ
 %����ֵ
    Theta_b2A1=0;
    Theta_b2A2=0;
    Theta_bA1=0;
    Theta_bA2=0;
    Ba2=0;
    for i =1:N
        Theta_bA1=Theta_bA1+PA1(1,i)*BActivity(i);
        Theta_bA2=Theta_bA2+PA2(1,i)*BActivity(i);
        Theta_b2A1=Theta_b2A1+PA1(1,i)*(BActivity(i)^2);
        Theta_b2A2=Theta_b2A2+PA2(1,i)*(BActivity(i)^2);
        Ba2=Ba2+BActivity(i)^2;
    end
    Theta_bA1=Theta_bA1/N;
    Theta_bA2=Theta_bA2/N;
    Theta_b2A2=Theta_b2A2/N;
    Theta_b2A1=Theta_b2A1/N;
    Ba=sum(BActivity)/N;
    Ba2=Ba2/N;
    x=(1-(1-gama1)* rhoA1-(1-gama2)* rhoA2)*(Ba2-(1-gama1)*Theta_b2A1-(1-gama2)*Theta_b2A2);
    EH=sqrt(x)+Ba-(1-gama1)*Theta_bA1-(1-gama2)*Theta_bA2;
%      zan=mu/(mb*EH);%��ֵ�ı��ʽ������ֻ���ڦ´���0.001ʱ�������ϸ�ֵ���������ֵ
     bata_c(xun)= mu/(mb*EH);  
end
%�洢��η����ֵ
beta_MMCA(rep,:)=bata_c;
beta_ave=mean(beta_MMCA,1);
end
xzhou=((1:ter)+13)/7;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
plot(xzhou,beta_ave)%','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
% plot(xzhou,beta_num(2,:)','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
% plot(xzhou,beta_num(3,:),'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
% set(gcf,'DefaultTextInterpreter','latex');
% xlabel('$\gamma_v$','FontSize',15);ylabel('$\beta_c$','FontSize',15);   % 鍧愭爣杞磋В閲�?
% h=legend('$\gamma_1=1.2,\gamma_2=1$','$\gamma_1=1,\gamma_2=1$','$\gamma_1=1,\gamma_2=0.8$');
% set(h,'Interpreter','latex','FontSize',15)%,
% % save('beta_num.mat','beta_num')
% % saveas(gcf, 'save.fig')
% 
% 
% 
% 
