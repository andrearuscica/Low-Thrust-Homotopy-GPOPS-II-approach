function [dae] = CartDaeScaled(sol,setup)

t     = sol.time;
xx    = sol.state;
x     = xx(:,1);
y     = xx(:,2);
z     = xx(:,3);
vx    = xx(:,4);
vy    = xx(:,5);
vz    = xx(:,6);
m     = xx(:,7);

T     = sol.control(:,1);

alpha = sol.control(:,2);
beta  = sol.control(:,3);
Isp = 3000;


r  = sqrt( x.^2 + y.^2 + z.^2 );


% Thrust Acceleration

cos_beta = cos(beta);
one_r3 =  1./(r.^3);

a = (T./m); %check that in the main it is already divided by 1000.
ax = a.*cos(alpha).*cos_beta;
ay = a.*sin(alpha).*cos_beta;
az = a .* sin(beta);

xdot  = vx;
ydot  = vy;
zdot  = vz;
vxdot = -one_r3 .*x + ax;
vydot = -one_r3 .*y + ay;
vzdot = -one_r3 .*z + az;
mdot = -T./(Isp*9.81).*1000.*setup.vc;

dae = [xdot, ydot, zdot, vxdot, vydot, vzdot, mdot];