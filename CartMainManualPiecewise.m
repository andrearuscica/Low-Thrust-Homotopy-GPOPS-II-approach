clear all; 

%% Problem Setup
data = GetProblemData(1); % 1: Johnson, 2: Barbee, 3: Moscow, 4: Moscow original
pars = load("1.GTOC_Data/pars.mat").pars; % spacecraft parameters

% Configuration
nphases = 6;
starting_arc = 10;
piece = 'string'; %Insert the previous solution file name

% Load previous solution data
prev_solution = load(piece);
piece_prev_solution_state = prev_solution.solution(1).state(end,:);
last_piece_velocity = piece_prev_solution_state(4:6);
last_piece_mass = piece_prev_solution_state(7);

% Spacecraft parameters
auxdata.T = pars.SC.T;
auxdata.Isp = pars.SC.Isp;
auxdata.mu = pars.mu_sun;
mu = auxdata.mu;

% Scaling factors
setup.nphases = nphases;
setup.starting_arc = starting_arc;
setup.total_arcs = height(data.arcs);
setup.lc = 149597870.700; % Astronomical Unit in km
setup.mu = 132712440018;  % Standard gravitational parameter
setup.tc = sqrt(setup.lc^3/setup.mu);
setup.vc = setup.lc/setup.tc;
setup.ac = setup.lc/setup.tc^2;
setup.m_init = 1500;

% Arc data
auxdata.arcsTable = data.arcs;

% Initial mass
if starting_arc == 1
    w0 = 1.0e+03 * 1.4829 / setup.m_init;
else
    w0 = last_piece_mass;
end
auxdata.w0 = w0;

%% Define Bounds
% Time bounds
t0 = 0; tmin = t0;
% Mass bounds
wmin = 500/setup.m_init;

% State bounds
xmin  = -3; xmax  = 3;
ymin  = -3; ymax  = 3;
zmin  = -1; zmax  = 1;
vxmin = -3; vxmax = 3;
vymin = -3; vymax = 3;
vzmin = -1; vzmax = 1;

% Control bounds
Tmax = 0.135/1000/(setup.m_init*setup.ac); Tmin = 0;
alphamin = -pi; alphamax = pi; 
betamin = -pi/2; betamax = pi/2;

% Get initial position and velocity based on starting arc
if starting_arc == 0
    arc_ri = data.Launch_state_ToF(1:3)/setup.lc;
    arc_prev_vf = data.Launch_state_ToF(4:6)/setup.vc;
elseif starting_arc == 1
    arc_ri = table2array(auxdata.arcsTable(starting_arc,2))/setup.lc;
    arc_prev_vf = data.EarthArc/setup.vc;
else
    arc_ri = table2array(auxdata.arcsTable(starting_arc,2))/setup.lc;
    arc_prev_vf = last_piece_velocity;
end

%% Configure Phase Limits
for iphase = 1:nphases
    % Set time of flight
    if starting_arc == 0 && iphase == 1
        arc_rf = table2array(auxdata.arcsTable(iphase,2))/setup.lc;
        arc_tof = data.Launch_state_ToF(7)./setup.tc;
    else
        arc_rf = table2array(auxdata.arcsTable(iphase+starting_arc-1,3))/setup.lc;
        arc_tof = table2array(auxdata.arcsTable(iphase+starting_arc-1,6))/setup.tc;
    end

    tf = arc_tof; tfmin = tf; tmax = tf;
    initial_state = [arc_ri arc_prev_vf w0]';
    
    % State bounds
    during_state_min = [xmin ymin zmin vxmin vymin vzmin wmin]';
    during_state_max = [xmax ymax zmax vxmax vymax vzmax initial_state(end)]';
    
    % Set up phase limits
    num_nodes = 10;
    limits(iphase).nodes = num_nodes;
    limits(iphase).time.min = [0 tfmin];
    limits(iphase).time.max = [0 tmax];

    % Configure state constraints
    if iphase == 1
        limits(iphase).state.min(:,1) = initial_state;
        limits(iphase).state.max(:,1) = initial_state;
    else
        limits(iphase).state.min(:,1) = during_state_min;
        limits(iphase).state.max(:,1) = during_state_max;
    end
    
    limits(iphase).state.min(:,2) = during_state_min;
    limits(iphase).state.min(:,3) = during_state_min;
    limits(iphase).state.max(:,2) = during_state_max;
    limits(iphase).state.max(:,3) = during_state_max;
    
    % Control constraints
    limits(iphase).control.min = [Tmin alphamin betamin]';
    limits(iphase).control.max = [Tmax alphamax betamax]';

    limits(iphase).parameter.min = [];
    limits(iphase).parameter.max = [];
    limits(iphase).path.min = [];
    limits(iphase).path.max = [];
    
    % Event constraints
    limits(iphase).event.min = [arc_rf]';
    limits(iphase).event.max = [arc_rf]';
end

%% Configure Linkages
for ipair = 1:nphases-1
    linkages(ipair).left.phase = ipair;
    linkages(ipair).right.phase = ipair + 1;
    linkages(ipair).min = zeros(7,1);
    linkages(ipair).max = zeros(7,1);
end
setup.limits = limits;

%% Setup Initial Guess
% Load Yago's guess
initialguess = load("Johnson_SF_10_21323_maxrep10.mat").sol_guess;
range = starting_arc:(starting_arc+nphases-1);

figure;
hold on;
for iphase = range
    % Scale state variables
    guess_position = initialguess.state(:,1:3,iphase)./setup.lc;
    guess_velocity = initialguess.state(:,4:6,iphase)./setup.vc;
    guess_mass = initialguess.mass(:,:,iphase)/setup.m_init;
    
    % Set time and state for current phase
    phase_idx = iphase-starting_arc+1;
    guess(phase_idx).time = (initialguess.times(:,:,iphase)/setup.tc);
    guess(phase_idx).state = [guess_position, guess_velocity, guess_mass];

    % Calculate control variables
    Ts = initialguess.thrust(:,:,iphase)./1000./(setup.m_init*setup.ac);
    Txs = Ts(:,1); Tys = Ts(:,2); Tzs = Ts(:,3);
    normTs = sqrt(Txs.^2+Tys.^2+Tzs.^2);
    beta = asin(Tzs./normTs);
    alpha = atan2(Tys,Txs);
    
    guess(phase_idx).control = [normTs, alpha, beta];
    guess(phase_idx).parameter = [];

    % Plot trajectory
    plot3(guess_position(:,1), guess_position(:,2), guess_position(:,3), 'r', 'LineWidth', 1.5);
end

grid on;
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
%% Full GPOPS Solution without optimization
prev_sol = load('Johnson FULL no opt.mat').solution;
for iphase = [1:prev_sol(nphases).phase]
guess(iphase).state = prev_sol(starting_arc+iphase-1).state;
guess(iphase).time = prev_sol(starting_arc+iphase-1).time';
guess(iphase).control = prev_sol(starting_arc+iphase-1).control;
guess(iphase).parameter = [];
hold on
plot3(guess(iphase).state(:,1),guess(iphase).state(:,2),guess(iphase).state(:,3),'r')
end

%% Piecewise with previous solution for every arc but the last, which uses GPOPS no opt
prev_sol = load(piece).solution;
for iphase = 1:(nphases-1)
    guess(iphase).state = prev_sol(iphase+1).state;
    guess(iphase).time = prev_sol(iphase+1).time';
    guess(iphase).control = prev_sol(iphase+1).control;
    guess(iphase).parameter = [];
    
    % Plot trajectory and final point
    plot3(guess(iphase).state(:,1), guess(iphase).state(:,2), guess(iphase).state(:,3), 'r', 'LineWidth', 1.5);
    x = guess(iphase).state(end,1); y = guess(iphase).state(end,2); z = guess(iphase).state(end,3);
    plot3(x, y, z, 'o', 'LineWidth', 1.5);
end

% Use GPOPS solution for last phase
gpops_noOpt = load('Johnson FULL no opt.mat').solution;
guess(nphases).state = gpops_noOpt(starting_arc+nphases-1).state;
guess(nphases).time = gpops_noOpt(starting_arc+nphases-1).time';
guess(nphase).control = gpops_noOpt(starting_arc+nphases-1).control;
guess(nphases).parameter = [];

plot3(guess(nphases).state(:,1), guess(nphases).state(:,2), guess(nphases).state(:,3), 'y', 'LineWidth', 1.5);
x = guess(nphases).state(end,1); y = guess(nphases).state(end,2); z = guess(nphases).state(end,3);
plot3(x, y, z, 'o', 'LineWidth', 1.5);

%% Piecewise with previous solution for every arc but the last, which uses Yago's
prev_sol = load(piece).solution;
for iphase = [1:prev_sol(nphases).phase-1]
guess(iphase).state = prev_sol(iphase+1).state;
guess(iphase).time = prev_sol(iphase+1).time';
guess(iphase).control = prev_sol(iphase+1).control;
guess(iphase).parameter = [];
end
yago_guess = load("Johnson_SF_10_21323_maxrep10.mat").sol_guess;

    guess_position = yago_guess.state(:,1:3,nphases+starting_arc-1)./setup.lc;
    guess_velocity = yago_guess.state(:,4:6,nphases+starting_arc-1)./setup.vc;
    guess_mass = yago_guess.mass(:,:,nphases+starting_arc-1)/setup.m_init;
    guess(nphases).time = (yago_guess.times(:,:,nphases+starting_arc-1)/setup.tc);
    guess(nphases).state = [guess_position, guess_velocity, guess_mass];
Ts = yago_guess.thrust(:,:,nphases+starting_arc-1)./1000./(setup.m_init*setup.ac);
Txs = Ts(:,1); Tys = Ts(:,2); Tzs = Ts(:,3);
normTs = sqrt(Txs.^2+Tys.^2+Tzs.^2);
beta = asin(Tzs./normTs);
alpha = atan2(Tys,Txs);
guess(nphases).control = [normTs,alpha,beta];
guess(nphases).parameter = [];
plot3(guess(nphases).state(:,1),guess(nphases).state(:,2),guess(nphases).state(:,3),'r')

%% Solver Configuration
setup.name = 'Trial Multiple Arcs';
setup.funcs.dae = 'CartDaeScaled';
setup.funcs.cost = 'CartCostMultipleScaled';
setup.funcs.event = 'CartEvent';
setup.funcs.link = 'CartLinks';
setup.auxdata = auxdata;
setup.guess = guess;
setup.linkages = linkages;
setup.derivatives = 'automatic';
setup.parallel = 'no';
setup.autoscale = 'off'; % Remember to turn on for pseudospectral method
setup.solver = 'ipopt';
setup.method = 'hermite-simpson';

%% Solve Problem
tic;
output = DMG(setup);
gpopsCPU = toc;
disp('END');
fprintf('Time to compile: %f \n', gpopsCPU);
solution = output.solution;

%% Save Solution
filepath = fullfile('C:', 'Users', 'andre', 'Desktop', 'Andrea', 'ISAE', 'ResearchProject', 'Codes', 'GPOPS_Mod', 'Johnson_piecewise', 'Trial 9', '9.sol_1_to_6_LAGRANGE');
save(filepath, 'solution');

%% Configure Plotting
% Set LaTeX as the interpreter for text elements
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

%% Plot Solution Trajectory
figure;
hold on;
for iphase = 1:nphases
    plot3(solution(iphase).state(:,1), solution(iphase).state(:,2), solution(iphase).state(:,3), 'LineWidth', 1.5);
end
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
grid on; view(3);

%% Prepare Continuous Time Vector
time_continuous = [];
end_times = [];
last_endtime = 0;

for i = 1:nphases
    time_continuous = [time_continuous, solution(i).time + last_endtime];
    last_endtime = solution(i).time(end) + last_endtime;
    end_times = [end_times, last_endtime];
end

%% Plot Mass History
masses = [];
for iphase = 1:nphases
    masses = [masses; solution(iphase).state(:,end)*setup.m_init];
end

figure;
hold on; grid on;
plot(time_continuous*setup.tc/86400, masses, 'LineWidth', 1.5);
xline(end_times*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1);
xlabel('Time (days)', 'FontSize', 14);
ylabel('Mass (kg)', 'FontSize', 14);
title('Mass History', 'FontSize', 16);
ax = gca;
ax.FontSize = 12;
ax.GridLineStyle = '--';

%% Plot Thrust History
thrusts = [];
for iphase = 1:nphases
    thrusts = [thrusts; solution(iphase).control(:,1)*1000*(setup.m_init*setup.ac)];
end

figure;
hold on; grid on;
plot(time_continuous*setup.tc/86400, thrusts, 'LineWidth', 1.5);
xline(end_times*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1);
yline(0.135, '--', 'Color', 'g', 'LineWidth', 1);
xlabel('Time (days)', 'FontSize', 14);
ylabel('Thrust (N)', 'FontSize', 14);
title('Thrust History', 'FontSize', 16);
ax = gca;
ax.FontSize = 12;
ax.GridLineStyle = '--';

%% Plot Steering History
alphas = [];
betas = [];
for iphase = 1:nphases
    alphas = [alphas; solution(iphase).control(:,2)];
    betas = [betas; solution(iphase).control(:,3)];
end

figure;
hold on; grid on;
plot(time_continuous*setup.tc/86400, alphas, 'LineWidth', 1.5);
plot(time_continuous*setup.tc/86400, betas, 'Color', '#77AC30', 'LineWidth', 1.5);
xline(end_times*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1);
xlabel('Time (days)', 'FontSize', 14);
ylabel('Angles (rad)', 'FontSize', 14);
title('Steering History', 'FontSize', 16);
legend({'\alpha', '\beta'}, 'Location', 'best', 'FontSize', 12);
ax = gca;
ax.FontSize = 12;
ax.GridLineStyle = '--';