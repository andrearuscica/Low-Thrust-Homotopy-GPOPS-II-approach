function event = Event(sol,setup)

%arcsTable = problem_data.arcs;
t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;


pf = xf(1);
ff = xf(2);
gf = xf(3);
hf = xf(4);
kf = xf(5);
Lf = xf(6);%*4*pi;

% ----------------------------------------------------------------------- %
% Pre-computing frequently used values
% ----------------------------------------------------------------------- %
cosLf = cos(Lf); sinLf = sin(Lf);

W = 1 + ff*cosLf+gf*sinLf; 
r = pf/W;
s2 = 1 + hf^2+kf^2;
alpha2 = hf^2 - kf^2;

rxf = r/s2 * (cosLf + alpha2 * cosLf + 2*hf*kf*sinLf);
ryf = r/s2 * (sinLf - alpha2 * sinLf + 2*hf*kf*cosLf);
rzf = 2*r/s2 * (hf*sinLf - kf*cosLf);


event = [rxf, ryf, rzf]';
end
