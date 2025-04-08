%CAMBIARE
function [Mayer,Lagrange]=CartCostMultipleScaled(sol,setup)


t  = sol.time;
xfin = sol.terminal.state;
xf = xfin(1);
yf = xfin(2);
zf = xfin(3);
vxf = xfin(4);
vyf = xfin(5);
vzf = xfin(6);
wf = xfin(7);


%Mayer = 0;
Mayer = -wf^2;

Lagrange = zeros(size(t));


 

end