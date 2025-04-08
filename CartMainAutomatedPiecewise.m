data = GetProblemData(1); % 1: Johnson, 2: Barbee, 3: Moscow, 4: Moscow original
total_arcs = height(data.arcs); setup.total_arcs = total_arcs;

pars = load("1.GTOC_Data/pars.mat").pars; %spacecraft parameters
%SCALING
setup.lc  = 149597870.700; %e03;
setup.mu  = 132712440018;  %e09;
setup.tc  = sqrt(setup.lc^3/setup.mu);
setup.vc  = setup.lc/setup.tc;
setup.ac  = setup.lc/setup.tc^2;
setup.m_init = 1500;
%Spacecraft
auxdata.T = pars.SC.T; T = pars.SC.T;
auxdata.Isp = pars.SC.Isp;
auxdata.mu = pars.mu_sun;
mu = auxdata.mu;
wmin = 500/setup.m_init;
Tmax = 0.135/1000/(setup.m_init*setup.ac); Tmin = 0; %% over 1000 to correctly adimentionalise [m/s2 to km/s2];
alphamin = -pi; alphamax = pi; betamin = -pi/2; betamax = pi/2;


%limits
xmin  = -3;
xmax  = 3;
ymin  = -3;
ymax  = 3;
zmin  = -1;
zmax  = 1;
vxmin = -3;
vxmax = 3;
vymin = -3;
vymax = 3;
vzmin = -1;
vzmax = 1;
%arcs
arcsTable = data.arcs; %i for rows, ri = 2,rf=3,vi = 4, vf= 5, tof = 6
auxdata.arcsTable = arcsTable;

nphases = 8; setup.nphases = nphases;
baseName1 = '8.sol_';
baseName2 = '_to_';
baseName3 = '.mat';

last_not_automated = 11;

num_nodes = 20; 

for i = last_not_automated:45-nphases

starting_arc = i+1;

piece = sprintf('%s%d%s%d%s',baseName1, i, baseName2, i+nphases-1, baseName3)
piece_prev_solution_state = load(piece).solution(1).state(end,:);
last_piece_velocity = piece_prev_solution_state(4:6);
last_piece_mass = piece_prev_solution_state(7);

w0 = last_piece_mass; auxdata.w0 = w0;
arc_prev_vf = last_piece_velocity;

t0 = 0; tmin = t0; 

arc_ri = table2array(arcsTable(starting_arc,2))/setup.lc;
for iphase = [1:nphases]
    arc_rf = table2array(arcsTable(iphase+starting_arc-1,3))/setup.lc;
    arc_tof = table2array(arcsTable(iphase+starting_arc-1,6))/setup.tc;

    tf = arc_tof; tfmin = tf; tmax = tf;
    initial_state = [arc_ri arc_prev_vf w0]';
    
    during_state_min = [xmin ymin zmin vxmin vymin vzmin wmin]';
    during_state_max = [xmax ymax zmax vxmax vymax vzmax initial_state(end)]';
    final_state_min = [xmin ymin zmin vxmin vymin vzmin wmin]';
    final_state_max = [xmax ymax zmax vxmax vymax vzmax initial_state(end)]';

    limits(iphase).time.min = [0 tfmin];
    limits(iphase).time.max = [0 tmax];

if iphase == 1
    limits(iphase).state.min(:,1) = initial_state;
    limits(iphase).state.min(:,2) = during_state_min;
    limits(iphase).state.min(:,3) = during_state_min;
    
    limits(iphase).state.max(:,1) = initial_state;
    limits(iphase).state.max(:,2) = during_state_max;
    limits(iphase).state.max(:,3) = during_state_max;
    else
    limits(iphase).state.min(:,1) = during_state_min;
    limits(iphase).state.min(:,2) = during_state_min;
    limits(iphase).state.min(:,3) = during_state_min;
        
    limits(iphase).state.max(:,1) = during_state_max;
    limits(iphase).state.max(:,2) = during_state_max;
    limits(iphase).state.max(:,3) = during_state_max;

end
    limits(iphase).nodes = num_nodes;
    limits(iphase).control.min = [Tmin alphamin betamin]';
    limits(iphase).control.max = [Tmax alphamax betamax]';

    limits(iphase).parameter.min = [];
    limits(iphase).parameter.max = [];
    
    limits(iphase).path.min = [];
    limits(iphase).path.max = [];
       
    limits(iphase).event.min = [arc_rf]';
    limits(iphase).event.max = [arc_rf]';
    
end
for ipair = 1:nphases-1
    %
    ipair = ipair;
    linkages(ipair).left.phase  = ipair;
    linkages(ipair).right.phase = ipair + 1;
    linkages(ipair).min = [0; 0; 0; 0; 0; 0; 0];
    linkages(ipair).max = [0; 0; 0; 0; 0; 0; 0];
    %
end
setup.limits = limits;

%guessPart

prev_sol = load(piece).solution;
for iphase = [1:prev_sol(nphases).phase-1]
guess(iphase).state = prev_sol(iphase+1).state;
guess(iphase).time = prev_sol(iphase+1).time';
guess(iphase).control = prev_sol(iphase+1).control;
guess(iphase).parameter = [];
hold on
plot3(guess(iphase).state(:,1),guess(iphase).state(:,2),guess(iphase).state(:,3),'r')
end

gpops_noOpt = load('Johnson FULL no opt.mat').solution;
guess(nphases).state = gpops_noOpt(starting_arc+nphases-1).state;
guess(nphases).time = gpops_noOpt(starting_arc+nphases-1).time';
guess(nphases).control = gpops_noOpt(starting_arc+nphases-1).control;
guess(nphases).parameter = [];
plot3(guess(nphases).state(:,1),guess(nphases).state(:,2),guess(nphases).state(:,3),'y')

%solverPart
setup.name = 'Trial Multiple Arcs';
setup.funcs.dae = 'CartDaeScaled';
setup.funcs.cost = 'CartCostMultipleScaled';
setup.funcs.event = 'CartEvent';
setup.auxdata = auxdata;
setup.guess = guess;
setup.funcs.link  = 'CartLinks';
setup.linkages = linkages;
setup.derivatives = 'automatic';
setup.parallel    = 'no';
setup.autoscale = 'off'; %Remember to on for pseudospectral method
setup.solver ='ipopt';
setup.method ='hermite-simpson'; %hermite-simpson

% ----------------------------------------------------------------------- %
% Solve problem and extract solution
% ----------------------------------------------------------------------- %
tic;
output = DMG(setup);
gpopsCPU = toc;
disp('END')
fprintf('Time to compile: %f \n', gpopsCPU)
solution = output.solution;

saveSolutionName = sprintf('%s%d%s%d%s',baseName1, i+1, baseName2, i+nphases, baseName3);
save(saveSolutionName, 'solution')

clear solution
clear guess



end