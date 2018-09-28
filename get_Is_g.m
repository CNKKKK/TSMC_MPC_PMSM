function [sys,x0,str,ts] = get_Is_g(t,x,u,flag)
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
sizes.NumOutputs     = 2;
sizes.NumInputs      = 9;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0.00001 0];     % 采样时间0.002s,500 Hz

% mdlOutputs

function sys=mdlOutputs(t,x,u)
global Is_A Is_B Is_C Us_A Us_B Us_C sita_u_g we Io_q
global Pi Po kesi delta wi Is_g Is_g_alphar Is_g_beta
global Us_max Ts Rs Rf Lf Cf flux frq

Us_A=u(1);Us_B=u(2);Us_C=u(3);
Is_A=u(4);Is_B=u(5);Is_C=u(6);
sita_u_g=u(7); we=u(8); Io_q=u(9);Ts=10e-6;Rs=0.958;
Rf=0.2;
Lf=0.002;
Cf=0.000021;

flux=0.1827;
frq=50;

Us_max=50;

wi=2*pi*frq;
Pi=1/2*Is_A*(1-2*wi^2*Cf*Lf)*(Us_A-Rf*Is_A)+1/2*Is_B*(1-2*wi^2*Cf*Lf)*(Us_B-Rf*Is_B)+1/2*Is_C*(1-2*wi^2*Cf*Lf)*(Us_C-Rf*Is_C);
Po=3/2*(flux*Io_q*we+Rs*Io_q^2);
kesi=Po/Pi;
delta=1/kesi*(-9*Rf*Rs*Io_q^2-9*flux*we*Rf)+9/4*Us_max^2;
Is_g=Us_max/(2*Rf)-sqrt(delta/(3*Rf));

 Is_g_alphar=Is_g*cos(sita_u_g);%电流电压同相位
 Is_g_beta=Is_g*sin(sita_u_g);
 
 sys(1)=Is_g_alphar;
 sys(2)=Is_g_beta;
