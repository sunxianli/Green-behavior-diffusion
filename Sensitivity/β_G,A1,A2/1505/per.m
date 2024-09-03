%�������
clc
clear
N=1000;
stp=200;
lamda1=0.6;
lamda2=0.6;
delta1=0.3;
delta2=0.3;
csmd_G=0.01*N;
csmd_A1=0.01*N;
csmd_A2=0.01*N;
gama1=1.5;
gama2=0.5;
mu=0.5;
ma=5;
mb=5;
%��ʼ����
termi=40;%�µ�ȡֵ����2
MC_rep=20;%�������
BETA_A1_NEW=zeros(MC_rep,termi);
BETA_A2_NEW=zeros(MC_rep,termi);
BETA_G_NEW=zeros(MC_rep,termi);
neighbora=1:N;%���������������Ĵ�ŵ�
neighborb=1:N;
BETA_A1_AVER =1:termi;%��ֵ����
BETA_A2_AVER =1:termi;
BETA_G_AVER =1:termi;
load('BActivity.mat')
load('AActivity.mat')
parfor rep = 1:MC_rep%�������
  %��ʼ��ͬ�Ħ�ֵѭ��
 BETA_A1=1:termi;%ÿ�η���Ľ��
 BETA_A2=1:termi;
 BETA_G=1:termi;%��ʼ����ѭ��
  for l = 1:termi
   disp(['��',num2str(rep),'�η���','��Ϊ',num2str(l)])
   beta_R=l/80;
%%�ȴ����²�Ľڵ㣬���²�ΪGʱ���ϲ�϶���A1,��������²�Ϊ0�����Զ���A2��U���Ͳ���Ҫ�����ˣ�Ȼ���ٲ����ϲ㣬ѡ1%��A2����
x=zeros(1,N);
m=zeros(1,N);
neighborA1=randi([1,N],1,2*csmd_G);
 for nA1 = 1:csmd_G;
          x(neighborA1(nA1))=1;%�²�G��Ӧ�ϲ�A1
          m(neighborA1(nA1))=1;
 end
 for nA2 =csmd_G+1:2*csmd_G;%%��ѡ���ϲ��A2���²��Զ���Ӧ0����S
          m(neighborA1(nA2))=2;
 end
   %��ʼ����ÿ���������з��棬����ʱ����stp=40
    for t = 1:stp
      n=1:N;
      n(:)=0;
      y=1:N;
      y(:)=0;
  
%m��x�ֱ�Ϊ���²��״̬����m��0��1��2��ɣ�x��0��1��ɣ�m��x��Ӧ����tʱ�̣�n��y��Ӧt+1ʱ�̵����²�
      A=zeros(N,N);%�ڽӾ���A�ϲ�
      B=zeros(N,N);%�²��ڽӾ���B
%�����ڽӾ���
    for i=1:N
    if rand(1)<AActivity(i)
        mnum=0;
        while(mnum<ma)
            randNum=floor(N*rand(1)+1);
            if A(i,randNum)==0 && i~=randNum
                    mnum=mnum+1;
                    A(i,randNum)=1;
                    A(randNum,i)=1;
            end
        end
    end
             if rand(1)<BActivity(i)
        mnum2=0;
        while(mnum2<mb)
            randNum=floor(N*rand(1)+1);
            if B(i,randNum)==0 && i~=randNum
                    mnum2=mnum2+1;
                   B(i,randNum)=1;
                    B(randNum,i)=1;
        end
    end
end
end

%%��ʼ�������²����ϴ���
      for i = 1:N
          %�ϲ㴫��
        if m(i)==0 ;%(U,ϵͳֻ�д���US������ʱx�ض�Ϊ0��US�ɱ�ΪA1S,A2S��
             for j = 1:N;
                 p1=rand(1,1);
                 p2=rand(1,1);  
            if (A(j,i)==1)&&(m(j)==2)&&(p2<lamda2);
                   n(i)=2;
                   break;
              elseif (A(j,i)==1)&&(m(j)==1)&&(p1<lamda1);
                   n(i)=1;
              end
            end 
        elseif m(i)==1;%A1��ϵͳ�д���A1S,A1G��A1S�ɱ�ΪUS��A1G�ɱ�UG 
            p1=rand(1,1);
                if p1<delta1;
                    n(i)=0;%US
                else
                    n(i)=1;
                end
        else %m(i)==2;%A2��ֻ����A2S,����ʱx��ȻΪ0��A2S�ɱ�ΪUS
            p2=rand(1,1);
            if p2<delta2;
              n(i)=0;  
            else
                n(i)=2;
            end
        end
         %�²㴫��
        if x(i)==0;%S״̬,ϵͳ�д���US,A1S,A2S��A1S����A1G��US�ȱ��UG,�ٱ��A1G;A2S�ȱ��A2G,�ٱ��A1G
          if n(i)==1;%A1S���ȴ����ӵ�A1S���ٴ�����ı��
               for j = 1:N
                   p1=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p1<gama1*beta_R);
                    y(i)=1;
                    break;
                end
               end
          elseif n(i)==0;%US���A1G
              for j = 1:N
                   p0=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p0<beta_R);
                    y(i)=1;
                    n(i)=1;
                        break;
                end
              end
          else n(i)==2;%A2S���ȱ��A2G,�ٱ��A1G
                for j = 1:N
                   p2=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p2<gama2*beta_R);
                    y(i)=1;
                    n(i)=1;
                        break;
                end
              end
          end
        else  x(i)==1;%G״̬��ϵͳ��ֻ����A1G,�����м�״̬UG,A1G�ɱ�ΪA1S,UG�ɱ�ΪUS��A1G
           if n(i)==1;%A1G
           p1=rand(1,1);
              if p1<mu;
                y(i)=0;   %A1G���ߣ����A1S
            else
                y(i)=1;    
              end
           else%if n(i)==0; %UG
                  p0=rand(1,1);
                  if p0<mu;%US
                      y(i)=0; 
                  else
                       y(i)=1; %A1G
                       n(i)=1;
                  end
           end
        end
      end
      %����m��x����
      m=n;
      x=y;  
    end
 %ͳ��ÿ���еĲ�ͬ�ڵ��������x����ֻ��0��1�����Ժܺü��㣻m����0��1��2Ҫ����һ�㣬�ȼ���m�еĽڵ���
 %for i=1:N
    A1=length(find(m(:)==1)); %����A1������
    A2=length(find(m(:)==2)); %����A2������             
    BETA_A1(l)=sum(A1)/N;%��l�η����У�A1���ܶ�
    BETA_A2(l)=sum(A2)/N;%��l�η����У�A2���ܶ�
    BETA_G(l)=sum(x)/N;%��l�η����У�G���ܶ�
  end
%����ÿ�η���Ľ��
  BETA_A1_NEW(rep,:)=BETA_A1;
  BETA_A2_NEW(rep,:)=BETA_A2;
  BETA_G_NEW(rep,:)=BETA_G;
end
%�����η����ƽ��ֵ
for ave = 1:termi
  BETA_A1_AVER(1,ave)=sum(BETA_A1_NEW(:,ave))./MC_rep;
  BETA_A2_AVER(1,ave)=sum(BETA_A2_NEW(:,ave))./MC_rep;
  BETA_G_AVER(1,ave)=sum(BETA_G_NEW(:,ave))./MC_rep;
end
xzhou=(1:termi)/80;
hold on;
box on;
grid off;
 set(gca,'Fontsize',15);
plot(xzhou,BETA_G_AVER','-o','color',[77/256 133/256 189/256]);
plot(xzhou,BETA_A1_AVER','-^','color',[247/256 144/256 61/256]);
plot(xzhou,BETA_A2_AVER,'-v','color',[89/256 169/256 90/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',15);ylabel('proportion','FontSize',15);   % ���������
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,
save('BETA_G_AVER.mat','BETA_G_AVER');
save('BETA_A1_AVER.mat','BETA_A1_AVER');
save('BETA_A2_AVER.mat','BETA_A2_AVER');
saveas(gcf, 'mc.fig');