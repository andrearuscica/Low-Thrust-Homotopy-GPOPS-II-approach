function event = CartEvent(sol,setup)

iphase = sol.phase;
t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xfin = sol.terminal.state;

xf = xfin(1);
yf = xfin(2);
zf = xfin(3);
vxf = xfin(4);
vyf = xfin(5);
vzf = xfin(6);

wf = xfin(7);


%starting_arc = setup.starting_arc;
%total_arcs = setup.total_arcs;

% if iphase + starting_arc - 1 == total_arcs %%Rendez-vous conditions
%     event = [xf, yf, zf, vxf, vyf, vzf]';
% else
    event = [xf, yf, zf]';
% end

end
