%CAMBIARE
function [Mayer,Lagrange]=Cost(sol,setup)

t  = sol.time;
xf = sol.terminal.state;
pf = xf(1);
ff = xf(2);
gf = xf(3);
hf = xf(4);
kf = xf(5);
Lf = xf(6);
wf = xf(7);
mu = 1.3272e11;

q = 1+ff*cos(Lf)+gf*sin(Lf);
r = pf/q;
alpha2 = hf*hf-kf*kf;
chi = sqrt(hf*hf+kf*kf);
s2 = 1+chi*chi;

 Mayer = -wf;
 %Mayer = 0;
 Lagrange = zeros(size(t));

end