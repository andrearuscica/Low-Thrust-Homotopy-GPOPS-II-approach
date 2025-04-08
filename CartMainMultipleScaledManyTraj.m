clear all; 
data = GetProblemData(1); % 1: Johnson, 2: Barbee, 3: Moscow, 4: Moscow original
total_arcs = height(data.arcs); setup.total_arcs = total_arcs;

pars = load("1.GTOC_Data/pars.mat").pars; %spacecraft parameters

nphases = 5;
starting_arc = 6; setup.starting_arc = starting_arc;

%Spacecraft
auxdata.T = pars.SC.T; T = pars.SC.T;
auxdata.Isp = pars.SC.Isp;
auxdata.mu = pars.mu_sun;
mu = auxdata.mu;

%SCALING
setup.lc  = 149597870.700; %e03;
setup.mu  = 132712440018;  %e09;
setup.tc  = sqrt(setup.lc^3/setup.mu);
setup.vc  = setup.lc/setup.tc;
setup.ac  = setup.lc/setup.tc^2;
setup.m_init = 1500;

%arcs
arcsTable = data.arcs; %i for rows, ri = 2,rf=3,vi = 4, vf= 5, tof = 6
auxdata.arcsTable = arcsTable;
if starting_arc == 0
    w0 = 1;
else
    w0 = table2array(arcsTable(starting_arc,7))/setup.m_init;
end
auxdata.w0 = w0;

t0 = 0; tmin = t0; 
wmin = 500/setup.m_init;

% rmin  = 0.5*AU;
% rmax  = 3*AU;
xmin  = -5;
xmax  = 5;
ymin  = -5;
ymax  = 5;
zmin  = -2;
zmax  = 2;
vxmin = -5;
vxmax = 5;
vymin = -5;
vymax = 5;
vzmin = -1;
vzmax = 1;



Tmax = 0.135/1000/(setup.m_init*setup.ac); Tmin = 0; % over 1000 to correctly adimentionalise [m/s2 to km/s2];
alphamin = -pi; alphamax = pi; betamin = -pi/2; betamax = pi/2;


if starting_arc == 0
    arc_ri = data.Launch_state_ToF(1:3)/setup.lc;
    arc_prev_vf =  data.Launch_state_ToF(4:6)/setup.vc;
    
elseif starting_arc == 1
        arc_ri = table2array(arcsTable(starting_arc,2))/setup.lc;
        arc_prev_vf = data.EarthArc/setup.vc;
else
        arc_ri = table2array(arcsTable(starting_arc,2))/setup.lc;
        arc_prev_vf = table2array(arcsTable(starting_arc-1,5))/setup.vc;
end

for iphase = [1:nphases]
    if starting_arc == 0 && iphase == 1
    arc_rf = table2array(arcsTable(iphase,2))/setup.lc;
    arc_tof = data.Launch_state_ToF(7)./setup.tc;
    else
    arc_rf = table2array(arcsTable(iphase+starting_arc-1,3))/setup.lc;
    arc_tof = table2array(arcsTable(iphase+starting_arc-1,6))/setup.tc;
    end


    
    

    tf = arc_tof; tfmin = tf; tmax = tf;
     initial_state = [arc_ri arc_prev_vf w0]';
    
   
    during_state_min = [xmin ymin zmin vxmin vymin vzmin wmin]';
    during_state_max = [xmax ymax zmax vxmax vymax vzmax initial_state(end)]';
    final_state_min = [xmin ymin zmin vxmin vymin vzmin wmin]';
    final_state_max = [xmax ymax zmax vxmax vymax vzmax initial_state(end)]';

    num_nodes = 20; limits(iphase).nodes = num_nodes;
    
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
    
    limits(iphase).control.min = [Tmin alphamin betamin]';
    limits(iphase).control.max = [Tmax alphamax betamax]';

    limits(iphase).parameter.min = [];
    limits(iphase).parameter.max = [];
    
    limits(iphase).path.min = [];
    limits(iphase).path.max = [];
       
    % event_tol = (1000/setup.lc)^2;
    % limits(iphase).event.min = [-sqrt(event_tol),-sqrt(event_tol),-sqrt(event_tol) , event_tol]';
    % limits(iphase).event.max = [sqrt(event_tol), sqrt(event_tol),sqrt(event_tol) , event_tol]';
    limits(iphase).event.min = [arc_rf]';
    limits(iphase).event.max = [arc_rf]';
    % 
    % if iphase + starting_arc - 1 == total_arcs
    %       limits(iphase).event.min = [data.Asteroids_CSV(end,1:3)./setup.lc, data.Asteroids_CSV(end,4:6)./setup.vc]';
    %       limits(iphase).event.max = [data.Asteroids_CSV(end,1:3)./setup.lc, data.Asteroids_CSV(end,4:6)./setup.vc]';
    % end
    % 
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

%% Yago's guess
initialguess = load("Johnson_SF_10_21323_maxrep10.mat").sol_guess;
range = [starting_arc:starting_arc+nphases-1];
hold on
for iphase = range
    guess_position = initialguess.state(:,1:3,iphase)./setup.lc;
    guess_velocity = initialguess.state(:,4:6,iphase)./setup.vc;
    guess_mass = initialguess.mass(:,:,iphase)/setup.m_init;
    guess(iphase-starting_arc+1).time = (initialguess.times(:,:,iphase)/setup.tc);
    guess(iphase-starting_arc+1).state = [guess_position, guess_velocity, guess_mass];

    Ts = initialguess.thrust(:,:,iphase)./1000./(setup.m_init*setup.ac);
    Txs = Ts(:,1); Tys = Ts(:,2); Tzs = Ts(:,3);
    normTs = sqrt(Txs.^2+Tys.^2+Tzs.^2);
    beta = asin(Tzs./normTs);
    alpha = atan2(Tys,Txs);
    guess(iphase-starting_arc+1).control = [normTs,alpha,beta];
    guess(iphase-starting_arc+1).parameter = [];

    %plotting
    plot3(guess_position(:,1),guess_position(:,2),guess_position(:,3),'r','LineWidth',1.5)
    hold on
end


hold on
grid on

xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]')

%% previous solution as initial guess
prev_sol = load('Johnson FULL no opt.mat').solution;
for iphase = [1:prev_sol(end).phase]
guess(iphase).state = prev_sol(iphase).state;
guess(iphase).time = prev_sol(iphase).time';
guess(iphase).control = prev_sol(iphase).control;
guess(iphase).parameter = [];
hold on
plot3(guess(iphase).state(:,1),guess(iphase).state(:,2),guess(iphase).state(:,3),'r')
end
%% Solver
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
%C:\Users\andre\Desktop\Andrea\ISAE\ResearchProject\Codes\GPOPS_Mod

%% LATEX INTERPRETER
% Set LaTeX as the interpreter for text elements
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

%% PLOT SOLUTION
%figure()
for iphase = [1:nphases]
   hold on
   plot3(solution(iphase).state(:,1),solution(iphase).state(:,2),solution(iphase).state(:,3),'LineWidth',1.5)
    alphas = solution(iphase).control(:,2);
    betas = solution(iphase).control(:,3);
    Ts = solution(iphase).control(:,1);
    Txs = Ts.*cos(alphas).*cos(betas); Tys = Ts.*sin(alphas).*cos(betas); Tzs = Ts.*sin(betas);
   %quiver3(solution(iphase).state(:,1),solution(iphase).state(:,2),solution(iphase).state(:,3),Txs,Tys,Tzs);

end
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]')
grid on; view(3);

%% TIME CONTINUITY
last_endtime = 0;
time_continuous = [];
end_times = [];
for i = [1:nphases]
time_continuous = [time_continuous, solution(i).time + last_endtime];
last_endtime = solution(i).time(end)+last_endtime;
end_times =[end_times, last_endtime];
end

%% PLOT MASS HISTORY
masses = [];
for iphase = 1:nphases
    masses = [masses; solution(iphase).state(:,end)*1500];
end
figure(); hold on; grid on
plot(time_continuous*setup.tc/86400,masses,'LineWidth',1.5)
xline(end_times*setup.tc/86400,'--', 'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex')
xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Mass$ (kg)', 'FontSize', 14, 'Interpreter', 'latex');
title('Mass History', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';
%saveas(gcf, 'full_noOPT_YagoGuess20_Mass_History.png');

%% PLOT THRUST HISTORY
thrusts = [];
for iphase = [1:nphases]
    thrusts = [thrusts; solution(iphase).control(:,1)*1000*(setup.m_init*setup.ac)];
end
figure(); hold on; grid on
plot(time_continuous*setup.tc/86400,thrusts,'LineWidth',1.5)
xline(end_times*setup.tc/86400,'--' ,'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex')

xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Thrust$ (N)', 'FontSize', 14, 'Interpreter', 'latex');
%yline(mean(thrusts),'-','Color','g','LineWidth',1.5,'Interpreter', 'latex')
title('Thrust History', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';
%saveas(gcf, 'full_02_noOPT_YagoGuess20_Thrust_History_WithAVG.png');

%% PLOT STEERING HISTORY
alphas = [];
betas = [];
for iphase = [1:nphases]
    alphas = [alphas; solution(iphase).control(:,2)];
    betas = [betas; solution(iphase).control(:,3)];
end


figure(); 
plot(time_continuous*setup.tc/86400,alphas,'LineWidth',1.5)
hold on; grid on;
plot(time_continuous*setup.tc/86400,betas,'Color','#77AC30','LineWidth',1.5)
xline(end_times*setup.tc/86400,'--', 'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex')


xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Angles$ (rad)', 'FontSize', 14, 'Interpreter', 'latex');
title('Steering History', 'FontSize', 16, 'Interpreter', 'latex');
legend({'$\alpha$','$\beta$'}, 'Location', 'best', 'FontSize', 12);
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';
saveas(gcf, 'full_02_noOPT_YagoGuess20_Steering_History.png');
