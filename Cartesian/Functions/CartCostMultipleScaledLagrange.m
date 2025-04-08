function [Mayer,Lagrange]=CartCostMultipleScaledLagrange(sol,setup)

t  = sol.time;
xfin = sol.terminal.state;
xf = xfin(1);
yf = xfin(2);
zf = xfin(3);
vxf = xfin(4);
vyf = xfin(5);
vzf = xfin(6);
wf = xfin(7);

u = sol.control;
T = u(:,1);
%Mayer = -wf^2;
Mayer = 0;



Lagrange = zeros(size(t));
weight = 0.1;
Lagrange = Lagrange + weight.*T.^2;

 

end