%-------------------------------------
% BEGIN: function Dae.m
%-------------------------------------
function dae = Dae(sol,setup);

auxdata = setup.auxdata;
Isp = auxdata.Isp; % [s]
mu = auxdata.mu; 

t = sol.time;
p = sol.state(:,1); 
f = sol.state(:,2); 
g = sol.state(:,3); 
h = sol.state(:,4);
k = sol.state(:,5);
L = sol.state(:,6);
w = sol.state(:,7);

u = sol.control;
ur = u(:,2);
ut = u(:,3);
un = u(:,4);

T = u(:,1);  %MODIFY THE THRUST DEFINITION


% ----------------------------------------------------------------------- %
% Thrust acceleration
% ----------------------------------------------------------------------- %
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;


DeltaT1 = (T./(1000*w)).*ur;
DeltaT2 = (T./(1000*w)).*ut;
DeltaT3 = (T./(1000*w)).*un;

% ----------------------------------------------------------------------- %
% Total acceleration
% ----------------------------------------------------------------------- %
Delta1 = DeltaT1;
Delta2 = DeltaT2;
Delta3 = DeltaT3;

% ----------------------------------------------------------------------- %
% Differential equations of motion
% ----------------------------------------------------------------------- %
dp = (2*p./q).*sqrt(p./mu).*Delta2;

df =  sqrt(p./mu).*sin(L).*Delta1 ...
     +sqrt(p./mu).*(1./q).*((q+1).*cos(L)+f).*Delta2 ...
     -sqrt(p./mu).*(g./q).*(h.*sin(L)-k.*cos(L)).*Delta3;

dg = -sqrt(p./mu).*cos(L).*Delta1 ...
     +sqrt(p./mu).*(1./q).*((q+1).*sin(L)+g).*Delta2 ...
     +sqrt(p./mu).*(f./q).*(h.*sin(L)-k.*cos(L)).*Delta3;

dh = sqrt(p./mu).*(s2.*cos(L)./(2*q)).*Delta3;

dk = sqrt(p./mu).*(s2.*sin(L)./(2*q)).*Delta3;

dL = sqrt(p./mu).*(1./q).*(h.*sin(L)-k.*cos(L)).*Delta3...
     +sqrt(mu.*p).*((q./p).^2);

dw = -T./(Isp*9.81);

path = (ur.^2 + ut.^2 + un.^2) - 1; 
dae = [dp,df,dg,dh,dk,dL,dw, path];



