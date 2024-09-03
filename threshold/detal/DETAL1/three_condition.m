%时变的A1A2，数值解不需要多次仿真
clc
clear
%参数设置
N=1000;
lamda1=0.6;
% lamda1=[0.1,0.25,0.4];
% mun=length(lamda1);
lamda2=0.6;
delta1=[0.1,0.7];
delta2=0.5
mun=length(delta1);
gama1=1.2;%衰减系数，当其为0时，知道信息的个体是完全免疫的。在这里若等于1，
gama2=0.8;
ma=5;%AD网络每次连接的边数
mb=5;
ter=20;% mu的取值个数
beta_MMCA=zeros(mun,ter);
stp=100;%时间长度
%定义矩阵
load AActivity;
load BActivity;
for l1=1:mun
  %  delta2=delta1(l1);
for xun = 1:ter
   mu=xun/ter;
   disp(['mu=' num2str(mu) ])
   
            PA2S=0.01*ones(1,N);
            PA1S=0.01*ones(1,N);

            PA2S_UPDATE=zeros(1,N);
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);
            rA2=zeros(1,N);

            RA1=zeros(N,N);
            RA2=zeros(N,N);

            %生成邻接矩阵
         
%计算马氏链
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
                    
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1(l1))+(1-PA1S(1,i)-PA2S(1,i))*(rA2(1,i)-rA1(1,i)*rA2(1,i));  
 PA2S_UPDATE(1,i)=PA2S(1,i)*(1-delta2)+(1-PA1S(1,i)-PA2S(1,i))*(1-rA2(1,i)); 
                end
                PA1S= PA1S_UPDATE;
                PA2S=PA2S_UPDATE;
            end       
%求阈值的图像时，当系统到达稳态时，G的值为0.所以A1=A1S
  PA1=PA1S;
  PA2=PA2S;
    
 rhoA1=sum(PA1)/N;
 rhoA2=sum(PA2)/N;%求出在某个β下的均值
 %求阈值
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
    zan=mu/(mb*EH);%阈值的表达式，但是只有在β大于0.001时，他的上个值才是真的阈值
    bata_c(l1,xun)=zan ;  
end
end
hold on;
 box on;
 grid off;
set(gca,'Fontsize',15);
xzhou=(1:ter)/ter;
plot(xzhou,bata_c(1,:)','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,bata_c(2,:)','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
% plot(xzhou,bata_c(3,:),'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\mu$','FontSize',15);ylabel('$\beta_c$','FontSize',15);   % 坐标轴解释
h=legend('$\delta_1$=0.1,$\delta_2$=0.5','$\delta_1$=0.7,$\delta_2$=0.5');
set(h,'Interpreter','latex','FontSize',15)%,
 save('bata_c.mat','bata_c')
 saveas(gcf, 'bata_c.fig');






