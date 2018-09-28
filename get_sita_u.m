function [sys,x0,str,ts] = get_sita_u(t,x,u,flag)
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
sizes.NumOutputs     = 1;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0.00001 0];     % 采样时间0.002s,500 Hz

% mdlOutputs

function sys=mdlOutputs(t,x,u)
global Id_g Iq_g sita_e_g
global sita_u_g
global faid_g faiq_g Ld Lq  L flux

Id_g=u(2);Iq_g=u(1); sita_e_g=u(3);
Ld=L;Lq=L; L=0.000395;flux=0.1827;
faid_g=flux+Ld*Id_g;
faiq_g=Lq*Iq_g;

sita_u_g=sita_e_g+atan(faiq_g/faid_g)+pi/2;
sys(1)=sita_u_g;
