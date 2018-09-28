function [sys,x0,str,ts] = switch_state(t,x,u,flag)
switch flag,
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1,
    sys=[];
  case 2,
    sys=[];
  case 3,
    sys=mdlOutputs(t,x,u);
  case 4,
    sys=[];
  case 9,
    sys=[];
   otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
% mdlInitializeSizes
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 14;
sizes.NumInputs      = 15;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0.00001 0];     % ����ʱ��0.0001s,100000 Hz

% mdlOutputs

function sys=mdlOutputs(t,x,u)
 persistent buffer a b c d e f g h i j k l 

Us_A=u(1);Us_B=u(2); Us_C=u(3);

Io_A_g=u(13); Io_B_g=u(14); Io_C_g=u(15);


Ts=0.00001;Rf=0.5; Lf=0.0004; Cf=21e-6; Rs=0.958; Ls=0.000395;%%���޸�������ȡIs��K+1)��״̬����ϵ��Ҳ����Ӧ���޸�
K1=0;K2=0; %%%%%K1 reactive power K2 source current


Io_A_1=u(10); Io_B_1=u(11); Io_C_1=u(12);
Uc_A_1=u(7); Uc_B_1=u(8); Uc_C_1=u(9);%%%%%%kshik
Is_A_1=u(4); Is_B_1=u(5); Is_C_1=u(6);%%%%%%%%%%%%%%����Uc Is Io
%%%%%%ʱ�Ӳ���%%%%%%%%%%%%%%%%%%%%%��kʱ�̵Ŀ���״̬������k+1ʱ�̵�Vs,Is,Io,Uo,Ic,Uc
% if isempty(buffer)%%�ж���һʱ�̵Ŀ���״̬�Ƿ�Ϊ�գ�Ϊ��һ�β�������ƣ�
% a=1; b=0;c=0;d=0;e=1;f=0;
% g=1;h=1;i=1;j=0;k=0;l=0;
% end
% Vdc_prev=(a-d)*Uc_A+(b-e)*Uc_B+(c-f)*Uc_C;
% Uo_A=Vdc_prev*(g-j); Uo_B=Vdc_prev*(h-k); Uo_C=Vdc_prev*(1/3)*(i-l);%%ʹ��k-1ʱ�̵Ŀ���״̬��ȡkʱ�������ѹ
%  Uo_A_c=Uo_A-(Uo_A+Uo_B+Uo_C)/3;Uo_B_c=Uo_B-(Uo_A+Uo_B+Uo_C)/3;Uo_C_c=Uo_C-(Uo_A+Uo_B+Uo_C)/3;%%%%%%%%%%%%%������ģ��ѹ
%  Idc_prev=g*Io_A+h*Io_B+i*Io_C;%%ʹ��K-1ʱ�̵Ŀ���״̬��ȡkʱ�̵�����任������
%  Ic_A=(a-d)*Idc_prev;Ic_B=(b-e)*Idc_prev;Ic_C=(c-f)*Idc_prev; 
% 
% 
% Is_A_1=-0.0729*Uc_A+0.9188*Is_A+0.0729*Us_A+0.0729*Ic_A;%%%%Is(k+1)=Is(k)+Uc(k)+Us(k)+Ic(k);
% Is_B_1=-0.0729*Uc_B+0.9188*Is_B+0.0729*Us_B+0.0739*Ic_B;
% Is_C_1=-0.0729*Uc_C+0.9188*Is_C+0.0729*Us_C+0.0739*Ic_C;
%  Io_A_1=1/(Ls)*((Ls-Ts*Rs)*Io_A+Ts*Uo_A_c);
%  Io_B_1=1/(Ls)*((Ls-Ts*Rs)*Io_B+Ts*Uo_B_c);%%ʹ��k-1ʱ�̵Ŀ���״̬��ȡkʱ���������
%  Io_C_1=1/(Ls)*((Ls-Ts*Rs)*Io_C+Ts*Uo_C_c);
% Uc_A_1=0.9261*Uc_A+1.9431*Is_A+0.0739*Us_A-1.9505*Ic_A;
% Uc_B_1=0.9261*Uc_B+1.9431*Is_B+0.0739*Us_B-1.9505*Ic_B;%%Io(k+1)=Io(k)+Uo(k)
% Uc_C_1=0.9261*Uc_C+1.9431*Is_C+0.0739*Us_C-1.9505*Ic_C;
%%%%%ʱ�Ӳ���
Uc_alphar_1=(2/3)* (Uc_A_1-(1/2)*Uc_B_1-(1/2)*Uc_C_1);      Uc_beta_1=(2/3)*(sqrt(3)/2)*(Uc_B_1-Uc_C_1);
Is_alphar_1=(2/3)* (Is_A_1-(1/2)*Is_B_1-(1/2)*Is_C_1);      Is_beta_1=(2/3)*((sqrt(3)/2)*Is_B_1-(sqrt(3)/2)*Is_C_1);
Io_g_alphar=(2/3)* (Io_A_g-(1/2)*Io_B_g-(1/2)*Io_C_g);      Io_g_beta=(2/3)*((sqrt(3)/2)*Io_B_g-(sqrt(3)/2)*Io_C_g);

%judge sector%
Us_alphar=(2/3)* (Us_A-(1/2)*Us_B-(1/2)*Us_C);      Us_beta=(2/3)*(sqrt(3)/2)*(Us_B-Us_C);
    if  Us_beta>0     
        A0=1;
    else
        A0=0;
    end
    if  0.866*Us_alphar-0.5*Us_beta>0    
        A1=1;        % sin(sita+2*pi/3)
    else
        A1=0;
    end    
    if  -0.866*Us_alphar-0.5*Us_beta>0     
        A2=1;       % sin(sita-2*pi/3)
    else
        A2=0;
    end
    P1=4*A2+2*A1+A0;

    if  P1==1            
        sector=2;
    elseif  P1==2        
        sector=6;
    elseif  P1==3        
        sector=1;
    elseif  P1==4       
        sector=4;
    elseif  P1==5        
        sector=3;
    elseif  P1==6       
        sector=5;
    else    
    end
%%%%%%%%%%%%%%%%%********%%%%%%%%%%%end
%%%%%%%%%%%%%%%%%%%%%���������Is%%%%%%%%%%%%

% Pi=(1-8*pi^2*50^2*Cf*Lf)*(Is_A_1*(Us_A-Rf*Is_A_1)+Is_B_1*(Us_B-Rf*Is_A_1)+Is_C_1*(Us_C-Rf*Is_C_1));
% Po=3*Rs*5^2;
landa=0.9;
kesi=1-8*pi^2*60^2*Cf*Lf;
Is_g=((-kesi*220)+sqrt((kesi*220)^2-4*kesi*Rf*Rs*(30)^2/landa))/(-2*kesi*Rf);
%%������ѹͬ��λ
%  Is_g=2.4;
 Is_g_A=Is_g*sin(2*pi*60*t);
 Is_g_B=Is_g*sin(2*pi*60*t+2/3*pi);
 Is_g_C=Is_g*sin(2*pi*60*t-2/3*pi);
 Is_g_alphar=(2/3)* (Is_g_A-(1/2)*Is_g_B-(1/2)*Is_g_C);      Is_g_beta=(2/3)*((sqrt(3)/2)*Is_g_B-(sqrt(3)/2)*Is_g_C);

%  Pi_A=(1-8*pi*50^2*Cf*Lf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sector==1

    %������ѡ�� S1 S5 ����ѡ��S1 S2 S3 on
    Vdc=1*Uc_A_1+(-1)*Uc_B_1; 
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0; %%����kʱ�̵Ŀ���״̬��ȡIc��K+1) �Լ�Uo��k+1��
     Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
     Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);%%%%ʹ��kʱ�̵ĵ�����ѹ��ȡK+1ʱ�̵ĵ���
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    %��Us(k+1)����Us��k��
    QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q1=K1*QK+QZ+K2*QZ1;
     switch_stat1=100010111000;
    
    
         %������ѡ�� S1 S5 ����ѡ��S4 S5 S6 on\
         
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
       Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    %��Us(k+1)����Us��k��
    QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010000111;
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
     %������ѡ�� S1 S5 ����ѡ��S4 S2 S3 on
 
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S1 S5 ����ѡ��S4 S5 S3 on
   
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
   %������ѡ�� S1 S5 ����ѡ��S1 S5 S6 on
     
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S1 S5 ����ѡ��S1 S2 S6 on
      
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
          %������ѡ�� S1 S5 ����ѡ��S1 S3 S5 on
   
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010101010;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
           %������ѡ�� S1 S5 ����ѡ��S2 S4 S6 on
   
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          
          
              %������ѡ�� S1 S6 ����ѡ��S1 S2 S3 on
    Vdc=1*Uc_A_1+(-1)*Uc_C_1; 
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001111000;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S1 S6 ����ѡ��S4 S5 S6 on
 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001000111;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
     %������ѡ�� S1 S6 ����ѡ��S4 S2 S3 on
    
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S1 S6 ����ѡ��S4 S5 S3 on
 
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
   %������ѡ�� S1 S6 ����ѡ��S1 S5 S6 on
  
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S1 S6 ����ѡ��S1 S2 S6 on
  
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
          %������ѡ�� S1 S6 ����ѡ��S1 S3 S5 on
   
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
           %������ѡ�� S1 S6 ����ѡ��S2 S4 S6 on
         Vdc=1*Uc_A_1+(-1)*Uc_C_1; 
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
        
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
    
        %������ѡ�� S2 S6 ����ѡ��S4 S5 S6 on
    
      Vdc=1*Uc_B_1+(-1)*Uc_C_1;
  

    Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
          Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001000111;

          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S2 S6 ����ѡ��S1 S2 S3 on

     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    
    QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001111000;
     
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S2 S6 ����ѡ��S4 S2 S3 on

     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
   
    QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S2 S6 ����ѡ��S4 S5 S3 on

     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
   %������ѡ�� S2 S6 ����ѡ��S1 S5 S6 on

     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S2 S6 ����ѡ��S1 S2 S6 on

     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
          %������ѡ�� S2 S6 ����ѡ��S1 S3 S5 on

     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
           %������ѡ�� S2 S6 ����ѡ��S2 S4 S6 on

     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
   Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001010101;
     
          if(Q1<Q2)
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          
          a=floor(switch_stat1/100000000000);
          b=floor((switch_stat1-100000000000*a)/10000000000);
          c=floor((switch_stat1-100000000000*a-10000000000*b)/1000000000);
          d=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c)/100000000);
          e=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d)/10000000);
          f=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e)/1000000);
          g=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f)/100000);
          h=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g)/10000);
          i=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h)/1000);
          j=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i)/100);
          k=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j)/10);   
          l=floor(switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j-10*k);
          
           buffer=[a b c d e f g h i j k l Q1 Is_g_B];
           sys=buffer;
          

elseif sector==2

    %������ѡ�� S2 S4 ����ѡ��S1 S2 S3 on
    Vdc=1*Uc_B_1+(-1)*Uc_A_1;
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    %��Us(k+1)����Us��k��
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q1=K1*QK+QZ+K2*QZ1;
     switch_stat1=010100111000;
     
     %������ѡ�� S2 S4 ����ѡ��S4 S5 S6 on
    
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
              Idc=0;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100000111;
          if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
          end
          
     %������ѡ�� S2 S4 ����ѡ��S4 S2 S3 on
   
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S2 S4 ����ѡ��S4 S5 S3 on
   
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S2 S4 ����ѡ��S1 S5 S6 on
    
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S2 S4 ����ѡ��S1 S2 S6 on
      
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S2 S4 ����ѡ��S1 S3 S5 on
     
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;     
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S2 S4 ����ѡ��S2 S4 S6 on
  
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
      
           %������ѡ�� S2 S6 ����ѡ��S4 S5 S6 on
    Vdc=1*Uc_B_1+(-1)*Uc_C_1; 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001000111;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
              %������ѡ�� S2 S6 ����ѡ��S1 S2 S3 on
    
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001111000;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
     %������ѡ�� S2 S6 ����ѡ��S4 S2 S3 on
 
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S2 S6 ����ѡ��S4 S5 S3 on
  
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S2 S6 ����ѡ��S1 S5 S6 on
 
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S2 S6 ����ѡ��S1 S2 S6 on
 
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S2 S6 ����ѡ��S1 S3 S5 on
    
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S2 S6 ����ѡ��S2 S4 S6 on
   
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
          
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
  
              
          
                    %������ѡ�� S1 S6 ����ѡ��S4 S5 S6 on

    Vdc=1*Uc_A_1+(-1)*Uc_C_1; 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
          Idc=0;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001000111;
          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S1 S6 ����ѡ��S1 S2 S3 on
 
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001111000;
     
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S1 S6 ����ѡ��S4 S2 S3 on
  
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S1 S6 ����ѡ��S4 S5 S3 on
 
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S1 S6 ����ѡ��S1 S5 S6 on
 
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S1 S6 ����ѡ��S1 S2 S6 on
     
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S1 S6 ����ѡ��S1 S5 S3 on
       
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S1 S6 ����ѡ��S4 S2 S6 on
  
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
 
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          a=floor(switch_stat1/100000000000);
          b=floor((switch_stat1-100000000000*a)/10000000000);
          c=floor((switch_stat1-100000000000*a-10000000000*b)/1000000000);
          d=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c)/100000000);
          e=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d)/10000000);
          f=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e)/1000000);
          g=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f)/100000);
          h=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g)/10000);
          i=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h)/1000);
          j=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i)/100);
          k=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j)/10);   
          l=floor(switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j-10*k);
          
           
           buffer=[a b c d e f g h i j k l Q1 Is_g_B];
           sys=buffer;
           
elseif sector==3   

%������ѡ�� S2 S4 ����ѡ��S1 S2 S3 on
    Vdc=1*Uc_B_1+(-1)*Uc_A_1; 
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    %��Us(k+1)����Us��k��
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);    
    Q1=K1*QK+QZ+K2*QZ1;
     switch_stat1=010100111000;
     
         %������ѡ�� S2 S4 ����ѡ��S4 S5 S6 on
    
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
         Idc=0;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100000111;
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S2 S4 ����ѡ��S4 S2 S3 on
    
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S2 S4 ����ѡ��S4 S5 S3 on
   
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S2 S4 ����ѡ��S1 S5 S6 on
    
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S2 S4 ����ѡ��S1 S2 S6 on
 
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S2 S4 ����ѡ��S1 S5 S3 on
     
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S2 S4 ����ѡ��S4 S2 S6 on
    
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta); 
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
        
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          
          
     %������ѡ�� S2 S6 ����ѡ��S4 S5 S6 on
    Vdc=1*Uc_B_1+(-1)*Uc_C_1; 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc; 
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001000111;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end     
          
              %������ѡ�� S2 S6 ����ѡ��S1 S2 S3 on
    
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001111000;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
     %������ѡ�� S2 S6 ����ѡ��S4 S2 S3 on
     
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S2 S6 ����ѡ��S4 S5 S3 on
   
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S2 S6 ����ѡ��S1 S5 S6 on
 
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta); 
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S2 S6 ����ѡ��S1 S2 S6 on
     
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S2 S6 ����ѡ��S1 S5 S3 on
      
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S2 S6 ����ѡ��S4 S2 S6 on
         
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta); 
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010001010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
       
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻����
     %������ѡ�� S3 S4 ����ѡ��S4 S5 S6 on

    Vdc=1*Uc_C_1+(-1)*Uc_A_1; 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100000111;
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end    
     
          %������ѡ�� S3 S4 ����ѡ��S1 S2 S3 on
  
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100111000;
     
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S3 S4 ����ѡ��S4 S2 S3 on
 
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S3 S4 ����ѡ��S4 S5 S3 on
   
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S3 S4 ����ѡ��S1 S5 S6 on
  
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S3 S4 ����ѡ��S1 S2 S6 on
  
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S3 S4 ����ѡ��S1 S3 S5 on
       
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S3 S4 ����ѡ��S2 S4 S6 on
    
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
  
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          a=floor(switch_stat1/100000000000);
          b=floor((switch_stat1-100000000000*a)/10000000000);
          c=floor((switch_stat1-100000000000*a-10000000000*b)/1000000000);
          d=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c)/100000000);
          e=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d)/10000000);
          f=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e)/1000000);
          g=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f)/100000);
          h=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g)/10000);
          i=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h)/1000);
          j=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i)/100);
          k=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j)/10);   
          l=floor(switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j-10*k);
          
           buffer=[a b c d e f g h i j k l Q1 Is_g_B];
           sys=buffer;
elseif sector==4 
%������ѡ�� S3 S5 ����ѡ��S4 S5 S6 on
%������һ��������п��ܻ���ΪVdcС��0������ѭ������ȱ��Q1���ʽ�������������
    Vdc=1*Uc_C_1+(-1)*Uc_B_1;
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc; 
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q1=K1*QK+QZ+K2*QZ1;
     switch_stat1=001010000111;

%������ѡ�� S2 S4 ����ѡ��S1 S2 S3 on

    Vdc=1*Uc_B_1+(-1)*Uc_A_1; 
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    %��Us(k+1)����Us��k��
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100111000;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end   
     
     
     %������ѡ�� S2 S4 ����ѡ��S4 S5 S6 on
   
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
         Idc=0;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100000111;
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end     
     
     %������ѡ�� S2 S4 ����ѡ��S4 S2 S3 on
    
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S2 S4 ����ѡ��S4 S5 S3 on
    
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S2 S4 ����ѡ��S1 S5 S6 on
    
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S2 S4 ����ѡ��S1 S2 S6 on
  
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S2 S4 ����ѡ��S1 S3 S5 on
     
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S2 S4 ����ѡ��S2 S4 S6 on
     
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=010100010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end   
          
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
         
          
          
      
              %������ѡ�� S3 S5 ����ѡ��S1 S2 S3 on  ��%������ѡ�� S3 S5 ����ѡ��S4 S5 S6
              %on����������
    Vdc=1*Uc_C_1+(-1)*Uc_B_1; 
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010111000;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
     %������ѡ�� S3 S5 ����ѡ��S4 S2 S3 on
  
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S3 S5 ����ѡ��S4 S5 S3 on
 
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
   %������ѡ�� S3 S5 ����ѡ��S1 S5 S6 on

     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S3 S5 ����ѡ��S1 S2 S6 on

     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S3 S5 ����ѡ��S1 S3 S5 on
   
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S3 S5 ����ѡ��S2 S4 S6 on
  
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
        
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
    
          
          %������ѡ�� S3 S4 ����ѡ��S4 S5 S6 on
    Vdc=1*Uc_C_1+(-1)*Uc_A_1; 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc; 
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100000111;
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end   
     
          %������ѡ�� S3 S4 ����ѡ��S1 S2 S3 on
   
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100111000;
     
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S3 S4 ����ѡ��S4 S2 S3 on
  
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
           
      %������ѡ�� S3 S4 ����ѡ��S4 S5 S3 on
  
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
   %������ѡ�� S3 S4 ����ѡ��S1 S5 S6 on
  
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1; 
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S3 S4 ����ѡ��S1 S2 S6 on
     
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_B_1; 
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
          %������ѡ�� S3 S4 ����ѡ��S1 S3 S5 on
  
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S3 S4 ����ѡ��S2 S4 S6 on
        
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
        QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          a=floor(switch_stat1/100000000000);
          b=floor((switch_stat1-100000000000*a)/10000000000);
          c=floor((switch_stat1-100000000000*a-10000000000*b)/1000000000);
          d=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c)/100000000);
          e=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d)/10000000);
          f=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e)/1000000);
          g=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f)/100000);
          h=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g)/10000);
          i=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h)/1000);
          j=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i)/100);
          k=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j)/10);   
          l=floor(switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j-10*k);
          
           buffer=[a b c d e f g h i j k l Q1 Is_g_B];
           sys=buffer;
elseif sector==5 
     %������ѡ�� S3 S5 ����ѡ��S4 S5 S6 on  ����������
    Vdc=1*Uc_C_1+(-1)*Uc_B_1;
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=0;Ic_B_1=Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
      QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q1=K1*QK+QZ+K2*QZ1;
     switch_stat1=001010000111;
    

    %������ѡ�� S1 S5 ����ѡ��S1 S2 S3 on
   
    Vdc=1*Uc_A_1+(-1)*Uc_B_1;
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    %��Us(k+1)����Us��k��
         QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010111000;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end     
    
    %������ѡ�� S1 S5 ����ѡ��S4 S5 S6 on
 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc; 
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
         Idc=0;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
          QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010000111;
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end     
     
     %������ѡ�� S1 S5 ����ѡ��S4 S2 S3 on
  
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
          QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
     
      %������ѡ�� S1 S5 ����ѡ��S4 S5 S3 on
  
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
          QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
   %������ѡ�� S1 S5 ����ѡ��S1 S5 S6 on

     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010100011;
    
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S1 S5 ����ѡ��S1 S2 S6 on
   
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010110001;

          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S1 S5 ����ѡ��S1 S3 S5 on

     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010101010;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S1 S5 ����ѡ��S2 S4 S6 on
      
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_B_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          
             
              %������ѡ�� S3 S5 ����ѡ��S1 S2 S3 on
    Vdc=1*Uc_C_1+(-1)*Uc_B_1;  
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010111000;
    
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
     %������ѡ�� S3 S5 ����ѡ��S4 S2 S3 on

     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
      Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010011100;
    
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S3 S5 ����ѡ��S4 S5 S3 on

     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S3 S5 ����ѡ��S1 S5 S6 on
 
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010100011;
 
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S3 S5 ����ѡ��S1 S2 S6 on

     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010110001;
   
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S3 S5 ����ѡ��S1 S3 S5 on
   
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010101010;
   
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S3 S5 ����ѡ��S2 S4 S6 on
       
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
            
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������

 %������ѡ�� S3 S4 ����ѡ��S4 S5 S6 on
    Vdc=1*Uc_C_1+(-1)*Uc_A_1;  
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100000111;
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end     
     
          %������ѡ�� S3 S4 ����ѡ��S1 S2 S3 on
    
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100111000;
     
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S3 S4 ����ѡ��S4 S2 S3 on
    
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_B_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S3 S4 ����ѡ��S4 S5 S3 on
 
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S3 S4 ����ѡ��S1 S5 S6 on
    
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S3 S4 ����ѡ��S1 S2 S6 on

     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S3 S4 ����ѡ��S1 S3 S5 on
        
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_A_1+Io_C_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100101010;
  
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S3 S4 ����ѡ��S2 S4 S6 on
      
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=-Idc;Ic_B_1=0;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001100010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
            
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          a=floor(switch_stat1/100000000000);
          b=floor((switch_stat1-100000000000*a)/10000000000);
          c=floor((switch_stat1-100000000000*a-10000000000*b)/1000000000);
          d=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c)/100000000);
          e=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d)/10000000);
          f=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e)/1000000);
          g=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f)/100000);
          h=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g)/10000);
          i=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h)/1000);
          j=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i)/100);
          k=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j)/10);   
          l=floor(switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j-10*k);
          
           buffer=[a b c d e f g h i j k l Q1 Is_g_B];
           sys=buffer;
else

    
    %������ѡ�� S1 S5 ����ѡ��S1 S2 S3 on
    Vdc=1*Uc_A_1+(-1)*Uc_B_1;
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    %��Us(k+1)����Us��k��
    QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q1=K1*QK+QZ+K2*QZ1;
     switch_stat1=100010111000;
     
     
             %������ѡ�� S1 S5 ����ѡ��S4 S5 S6 on
    Vdc=1*Uc_A_1+(-1)*Uc_B_1;  
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
         Idc=0;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
    QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010000111;
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S1 S5 ����ѡ��S4 S2 S3 on
    
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S1 S5 ����ѡ��S4 S5 S3 on
   
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end

   %������ѡ�� S1 S5 ����ѡ��S1 S5 S6 on
   
     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1; 
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S1 S5 ����ѡ��S1 S2 S6 on
   
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S1 S5 ����ѡ��S1 S3 S5 on
        
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010101010;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S1 S5 ����ѡ��S2 S4 S6 on
     
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=Idc;Ic_B_1=-Idc;Ic_C_1=0;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100010010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
 %������ѡ�� S3 S5 ����ѡ��S4 S5 S6 on

    Vdc=1*Uc_C_1+(-1)*Uc_B_1; 
     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;     
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
               Idc=0;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010000111;
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
              %������ѡ�� S3 S5 ����ѡ��S1 S2 S3 on
  
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010111000;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
     %������ѡ�� S3 S5 ����ѡ��S4 S2 S3 on
 
     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
            
      %������ѡ�� S3 S5 ����ѡ��S4 S5 S3 on

     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          

   %������ѡ�� S3 S5 ����ѡ��S1 S5 S6 on

     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S3 S5 ����ѡ��S1 S2 S6 on
 
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S3 S5 ����ѡ��S1 S3 S5 on
        
     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S3 S5 ����ѡ��S2 S4 S6 on
        
     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=0;Ic_B_1=-Idc;Ic_C_1=Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=001010010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; 
              switch_stat1=switch_stat2;
          end
         
          
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
          
    
          
          %������ѡ�� S1 S6 ����ѡ��S1 S2 S3 on
    Vdc=1*Uc_A_1+(-1)*Uc_C_1;  
     Uo_A_1=1/3*(2-1-1)*Vdc;Uo_B_1=1/3*(-1+2-1)*Vdc;Uo_C_1=1/3*(-1-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001111000;
     
     if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
     end
     
     %������ѡ�� S1 S6 ����ѡ��S4 S2 S3 on

     Uo_A_1=1/3*(2*0-1-1)*Vdc;Uo_B_1=1/3*(-1*0+2-1)*Vdc;Uo_C_1=1/3*(-1*0-1+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001011100;
     
     if(Q1<Q2) 
     else
         Q1=Q2; switch_stat1=switch_stat2;
     end
             
      %������ѡ�� S1 S6 ����ѡ��S4 S5 S3 on
   
     Uo_A_1=1/3*(2*0-1*0-1)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001001110;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
    
   %������ѡ�� S1 S6 ����ѡ��S1 S5 S6 on

     Uo_A_1=1/3*(2-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001100011;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
     
           %������ѡ�� S1 S6 ����ѡ��S1 S2 S6 on
   
     Uo_A_1=1/3*(2-1-1*0)*Vdc;Uo_B_1=1/3*(-1+2-1*0)*Vdc;Uo_C_1=1/3*(-1-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_B_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001110001;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %������ѡ�� S1 S6 ����ѡ��S1 S3 S5 on

     Uo_A_1=1/3*(2-1*0-1)*Vdc;Uo_B_1=1/3*(-1+2*0-1)*Vdc;Uo_C_1=1/3*(-1-1*0+2)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_A_1+Io_C_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001101010;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
           %������ѡ�� S1 S6 ����ѡ��S2 S4 S6 on

     Uo_A_1=1/3*(2*0-1-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=Io_B_1;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);  
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001010101;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
          
          %����ѡ��S1S2S3������S4S5S6������Ч��һ�£��㲻������
                     %������ѡ�� S1 S6 ����ѡ��S4 S5 S6 on

     Uo_A_1=1/3*(2*0-1*0-1*0)*Vdc;Uo_B_1=1/3*(-1*0+2*0-1*0)*Vdc;Uo_C_1=1/3*(-1*0-1*0+2*0)*Vdc;
      Uo_A_1_c=Uo_A_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_B_1_c=Uo_B_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;Uo_C_1_c=Uo_C_1-(Uo_A_1+Uo_B_1+Uo_C_1)/3;%%������ģ��ѹ
     Idc=0;
     Ic_A_1=Idc;Ic_B_1=0;Ic_C_1=-Idc;
      Ic_alphar_1=(2/3)* (Ic_A_1-(1/2)*Ic_B_1-(1/2)*Ic_C_1);      Ic_beta_1=(2/3)*((sqrt(3)/2)*Ic_B_1-(sqrt(3)/2)*Ic_C_1);
    Is_alphar_2=-0.0248*Uc_alphar_1+0.9817*Is_alphar_1+0.0248*Us_alphar+0.0059*Ic_alphar_1;
     Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1_c);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1_c);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);
       Io_2_alphar=(2/3)* (Io_A_2-(1/2)*Io_B_2-(1/2)*Io_C_2);    Io_2_beta=(2/3)*((sqrt(3)/2)*Io_B_2-(sqrt(3)/2)*Io_C_2);
     QK=(abs(Is_alphar_2*Us_beta-Us_alphar*Is_beta_2));QZ=(abs(Io_g_alphar-Io_2_alphar)+abs(Io_g_beta-Io_2_beta));QZ1=abs(Is_alphar_2-Is_g_alphar)+abs(Is_beta_2-Is_g_beta);
    Q2=K1*QK+QZ+K2*QZ1;
     switch_stat2=100001000111;
     
          if(Q1<Q2) 
          else
              Q1=Q2; switch_stat1=switch_stat2;
          end
         
          a=floor(switch_stat1/100000000000);
          b=floor((switch_stat1-100000000000*a)/10000000000);
          c=floor((switch_stat1-100000000000*a-10000000000*b)/1000000000);
          d=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c)/100000000);
          e=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d)/10000000);
          f=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e)/1000000);
          g=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f)/100000);
          h=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g)/10000);
          i=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h)/1000);
          j=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i)/100);
          k=floor((switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j)/10);   
          l=floor(switch_stat1-100000000000*a-10000000000*b-1000000000*c-100000000*d-10000000*e-1000000*f-100000*g-10000*h-1000*i-100*j-10*k);
          
           buffer=[a b c d e f g h i j k l Q1 Is_g_B];
           sys=buffer;

end

        
       
   
     
    
