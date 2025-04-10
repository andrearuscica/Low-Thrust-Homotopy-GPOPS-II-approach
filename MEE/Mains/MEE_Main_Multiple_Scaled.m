clear all; clc
data = GetProblemData(1); % 1: Johnson, 2: Barbee, 3: Moscow, 4: Moscow original

pars = load("1.GTOC_Data/pars.mat").pars; %spacecraft parameters

nphases = 3;
starting_arc = 10;

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

fmin = -1; fmax = +1;
gmin = -1; gmax = +1;
hmin = -1; hmax = +1;
kmin = -1; kmax = +1;
%L is inside the for loop
wmin = 500/setup.m_init; wmax = w0;
urmin = -1; urmax = +1;
utmin = -1; utmax = +1;
uhmin = -1; uhmax = +1;
Tmin = 0; Tmax = T*1000000/1000/(setup.m_init*setup.ac); % over 1000 to correctly adimentionalise [m/s2 to km/s2];




% 
% if starting_arc == 0
%     arc_ri = data.Launch_state_ToF(1:3)/setup.lc;
%     arc_prev_vf =  data.Launch_state_ToF(4:6)/setup.vc;
% 
% elseif starting_arc == 1
%         arc_ri = table2array(arcsTable(starting_arc,2))/setup.lc;
%         arc_prev_vf = data.EarthArc/setup.vc;
% else
%         arc_ri = table2array(arcsTable(starting_arc,2))/setup.lc;
%         arc_prev_vf = table2array(arcsTable(starting_arc-1,5))/setup.vc;
% end

Lscale = 4*pi;

arc_ri = data.Launch_state_ToF(1:3);
arc_vi = data.Launch_state_ToF(4:6);
[a,e,i,RAAN,w,MA,nu] = orbitalElements(arc_ri, arc_vi, pars);
initialMEE = orbital2equinoctial([a,e,i,RAAN,w,nu]);
p0 = initialMEE(1); f0 = initialMEE(2); g0 = initialMEE(3); h0 = initialMEE(4); k0 = initialMEE(5); L0 = initialMEE(6)/Lscale;

initial_state = [p0 f0 g0 h0 k0 L0 w0]';

for iphase = [1:nphases]
    if starting_arc == 0 && iphase == 1
    arc_rf = table2array(arcsTable(iphase,2))/setup.lc;
    arc_tof = data.Launch_state_ToF(7)./setup.tc;
    else
    arc_rf = table2array(arcsTable(iphase+starting_arc-1,3))/setup.lc
    arc_tof = table2array(arcsTable(iphase+starting_arc-1,6))/setup.tc;
    
    Lmin = L0; Lmax = 9*2*pi;
    w0 = table2array(arcsTable(iphase+starting_arc-1,7))/setup.m_init;
    end


    
    

    tf = arc_tof; tfmin = tf; tmax = tf;
   
   
    during_state_min = [0.5*p0 fmin gmin hmin kmin Lmin-0.001 wmin]';
    during_state_max = [1.5*p0 fmax gmax hmax kmax Lmax initial_state(end)]'; 
    final_state_min = [0.5*p0 fmin gmin hmin kmin Lmin-0.001 wmin]';
    final_state_max = [1.5*p0 fmax gmax hmax kmax Lmax initial_state(end)]';

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
    
    limits(iphase).control.min = [Tmin urmin utmin uhmin]';
    limits(iphase).control.max = [Tmax urmax utmax uhmax]';

    limits(iphase).parameter.min = [];
    limits(iphase).parameter.max = [];
    
    limits(iphase).path.min = 0;
    limits(iphase).path.max = 0;
       
    % event_tol = (1000/setup.lc)^2;
    % limits(iphase).event.min = [-sqrt(event_tol),-sqrt(event_tol),-sqrt(event_tol) , event_tol]';
    % limits(iphase).event.max = [sqrt(event_tol), sqrt(event_tol),sqrt(event_tol) , event_tol]';
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

%% Yago's guess
initialguess = load("Johnson_10_21322_1h29min.mat").sol_guess;
%initialguess = load('Moscow_10_21322_Best_1h51min.mat').sol_guess;

range = [starting_arc:starting_arc+nphases-1];
hold on
for current_arc = range
    guess_state_OE = zeros(length(initialguess.state(:,:,current_arc)),6);
    guess_state_MEE = zeros(length(initialguess.state(:,:,current_arc)),7);
    for j = [1:length(initialguess.state(:,:,current_arc))]
        [a,e,i,RAAN,w,MA,nu] = orbitalElements(initialguess.state(j,1:3,current_arc), initialguess.state(j,4:6,current_arc), pars);
        guess_state_OE(j,:) = [a,e,i,RAAN,w,nu];
        guess_state_MEE(j,1:6) = orbital2equinoctial(guess_state_OE(j,:));
        guess_state_MEE(j,7) = initialguess.mass(j,:,current_arc)./setup.m_init;
        
    end
    guess_state_MEE(:,6) = guess_state_MEE(:,6);
    guess(current_arc-starting_arc+1).time = (initialguess.times(:,:,current_arc)/setup.tc);
    guess(current_arc-starting_arc+1).state = guess_state_MEE;


    r_vec = initialguess.state(:,1:3,current_arc);
    v_vec = initialguess.state(:,4:6,current_arc);
    
    normR = (r_vec(:,1).^2 + r_vec(:,2).^2 + r_vec(:,3).^2).^(1/2);
    normV = (v_vec(:,1).^2 + v_vec(:,2).^2 + v_vec(:,3).^2).^(1/2);
    
    ir = r_vec./normR;
    h_vec = cross(r_vec,v_vec,2);
    normH = (h_vec(:,1).^2 + h_vec(:,2).^2 + h_vec(:,3).^2).^(1/2);
    in = (h_vec)./normH;
    it = cross(in,ir,2);
    
    for i = [1:height(ir)]
        %each "page" of the matrix is a conversion matrix for a state
        Q_ECI2MEE(:,:,i) = [ir(i,:); it(i,:); in(i,:)];
    end
    
    t1 = initialguess.thrust(:,1,current_arc);
    t2 = initialguess.thrust(:,2,current_arc);
    t3 = initialguess.thrust(:,3,current_arc);
    magT = sqrt(t1.^2 + t2.^2 + t3.^2);
    T1 = initialguess.thrust(:,1)./magT;
    T2 = initialguess.thrust(:,2)./magT;
    T3 = initialguess.thrust(:,3)./magT;
    
    for i = [1:length(Q_ECI2MEE(1,1,:))]
        uMEE(:,:,i) = Q_ECI2MEE(:,:,i)*[T1(i) T2(i) T3(i)]';
    end
    
    %reverse to double check
    for i = [1:length(Q_ECI2MEE(1,1,:))]
        uECI(:,:,i) = Q_ECI2MEE(:,:,i)'*uMEE(:,:,i);
    end
    
    
    
    
    guess_control = zeros(length(initialguess.state(:,:,current_arc+1)),4);
    guess_control(:,1) = magT/1000/(setup.m_init*setup.ac);
    for i = 1:length(uMEE(1,1,:))
    guess_control(i,2:4) = uMEE(:,:,i)';
    end
    
    %initialguess.parameter(1:length(initialguess.phase.control(:,1)),1) = initialguess.parameter;
    
    guess(current_arc-starting_arc+1).time = initialguess.times(:,:,current_arc+1)./setup.tc;
    guess(current_arc-starting_arc+1).state = guess_state_MEE;
    guess(current_arc-starting_arc+1).control = guess_control;
    guess(current_arc-starting_arc+1).parameter = [];


   

    %plotting
    %plot3(guess_position(:,1),guess_position(:,2),guess_position(:,3),'r','LineWidth',1.5)
    hold on
    grid on

xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]')

p = guess_state_MEE(:,1);
f = guess_state_MEE(:,2);%*fscale;
g = guess_state_MEE(:,3);%*gscale;
h = guess_state_MEE(:,4);%*hscale;
k = guess_state_MEE(:,5);%*kscale;
L = guess_state_MEE(:,6);

q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;

% ----------------------------------------------------------------------- %
% Trajectory in Cartesian
% ----------------------------------------------------------------------- % 
X = (r./s2).*(cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L));
Y = (r./s2).*(sin(L)-alpha2.*sin(L)+2.*h.*k.*cos(L));
Z = (2.*r./s2).*(h.*sin(L)-k.*cos(L));
vX = -(1./s2).*sqrt(1./p).*(sin(L)+alpha2.*sin(L)-2.*h.*k.*cos(L)+g-2.*f.*h.*k+alpha2.*g);
vY = -(1./s2).*sqrt(1./p).*(-cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L)-f+2.*g.*h.*k+alpha2.*f);
vZ = (2./s2).*sqrt(1./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);

plot3(X,Y,Z)
end



%% previous solution
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
setup.funcs.dae = 'MEEDaeScaled4';
setup.funcs.cost = 'Cost';
setup.funcs.event = 'Event';
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
% At each phase, the times start from 0 so to plot  
% them in the same graph they need to be adjusted
times_continuous = zeros(num_nodes,1,nphases);
for i = [1:nphases]
    if i == 1
        times_continuous1(:,:,i) = solution(i).time';
    else
        times_continuous1(:,:,i) = solution(i).time' + times_continuous(end,end,i-1);
     end
end
%%
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

%% PLOT MASS HISTORY
figure()
for iphase = [1:nphases]
    hold on; grid on
    plot(times_continuous(:,:,iphase)*setup.tc/86400, solution(iphase).state(:,end)*1500,'LineWidth',1.5);
    xline(times_continuous(end,end,iphase)*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex');
    
end
xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Mass$ (kg)', 'FontSize', 14, 'Interpreter', 'latex');
title('Mass History', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';
saveas(gcf, 'Arc3_to_7Multiple30_60_last_two_VI_Mass_History.png');
%% PLOT THRUST HISTORY
figure()
for iphase = [1:nphases]
    hold on; grid on
    plot(times_continuous(:,:,iphase)*setup.tc/86400, solution(iphase).control(:,1)*1000*(setup.m_init*setup.ac),'LineWidth',1.5);
    xline(times_continuous(end,end,iphase)*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex');
    
end
xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Thrust$ (N)', 'FontSize', 14, 'Interpreter', 'latex');
title('Thrust History', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';
%saveas(gcf, 'Arc_27_to_30_08_05_Thrust_History.png');

%% ALPHA BETA ANGLE HISTORY
figure()
for iphase = [1:nphases]
    hold on; grid on
    plot(times_continuous(:,:,iphase)*setup.tc/86400, solution(iphase).control(:,2),'b','LineWidth',1.5);
    plot(times_continuous(:,:,iphase)*setup.tc/86400, solution(iphase).control(:,3),'r','LineWidth',1.5);
    xline(times_continuous(end,end,iphase)*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex');
    
end
xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Angles$ (rad)', 'FontSize', 14, 'Interpreter', 'latex');
title('Steering History', 'FontSize', 16, 'Interpreter', 'latex');
legend({'$\alpha$','$\beta$'}, 'Location', 'best', 'FontSize', 12);
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';
saveas(gcf, 'Arc3_to_7Multiple30_60_last_two_VI_Steering_History.png');

    %plotting
    plot3(guess_position(:,1),guess_position(:,2),guess_position(:,3),'r','LineWidth',1.5)
    hold on



hold on
%quiver3(initialguess.state(:,1), initialguess.state(:,2), initialguess.state(:,3),initialguess.thrust(:,1),initialguess.thrust(:,2),initialguess.thrust(:,3))
grid on

xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]')
