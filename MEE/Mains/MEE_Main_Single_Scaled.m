clear all; clc; %close all;

%global arc_index
arc_index = 13;
auxdata.arc_index = arc_index;

addpath("1.GTOC_Data")
addpath("0.Auxiliary_functions")
addpath("2.Initial_guesses")
pars = load("1.GTOC_Data/pars.mat").pars;
problem_data = load("0.Auxiliary_functions/problem_data.mat").problem_data;
AU = 1.496e8;
%Parameters
auxdata.T = pars.SC.T; T = pars.SC.T;
auxdata.Isp = pars.SC.Isp;
auxdata.mu = pars.mu_sun;
mu = auxdata.mu;
%
auxdata.mu = pars.mu_sun;
%SCALING
setup.lc  = pars.AU; %e03;
setup.mu  = 132712440018;  %e09;
setup.tc  = sqrt(setup.lc^3/setup.mu);
setup.vc  = setup.lc/setup.tc;
setup.ac  = setup.lc/setup.tc^2;
setup.m_init = 1500;

%arcs
arcsTable = problem_data.arcs; %i for rows, ri = 2,rf=3,vi = 4, vf= 5, tof = 6
arc_tof = table2array(arcsTable(arc_index,6))/setup.tc;
arc_prev_vf = table2array(arcsTable(arc_index-1,5));
w0 = table2array(arcsTable(arc_index,7))/setup.m_init;

auxdata.w0 = w0;
t0 = 0; tmin = t0; tf = arc_tof; tfmin = tf; tmax = tf;

MEE_Table = problem_data.MEE;

%fscale = 0.25; gscale = 0.05; hscale = 0.05; kscale = 0.03;
previous_arc_MEE = table2array(MEE_Table(arc_index-1,2));
p0 = previous_arc_MEE(1);
f0 = previous_arc_MEE(2);
g0 = previous_arc_MEE(3);
h0 = previous_arc_MEE(4);
k0 = previous_arc_MEE(5);
L0 = previous_arc_MEE(6)/(4*pi);

fmin = -1; fmax = +1;
gmin = -1; gmax = +1;
hmin = -1; hmax = +1;
kmin = -1; kmax = +1;
Lmin = L0; Lmax = 9;%*2*pi;
wmin = 500/setup.m_init; wmax = w0;
urmin = -1; urmax = +1;
utmin = -1; utmax = +1;
uhmin = -1; uhmax = +1;
Tmin = 0; Tmax = T/1000/(setup.m_init*setup.ac);

initial_state = [p0 f0 g0 h0 k0 L0 w0]';
%initial_state = load("prova_single_optimization_arc8_9.mat").prev_arc_sol';
during_state_min = [0.5*p0 fmin gmin hmin kmin Lmin-0.001 wmin]';
during_state_max = [1.5*p0 fmax gmax hmax kmax Lmax initial_state(end)]'; 
final_state_min = [0.5*p0 fmin gmin hmin kmin Lmin-0.001 wmin]';
final_state_max = [1.5*p0 fmax gmax hmax kmax Lmax initial_state(end)]';
% initial MEE dell'orbita del previous lambert;
% final a piacere;

r_event = table2array(arcsTable(arc_index,3))/setup.lc;

%% LIMITS
num_nodes = 60; limits(1).nodes = num_nodes;

limits(1).time.min = [0 tfmin];
limits(1).time.max = [0 tmax];

limits(1).state.min(:,1) = initial_state;
limits(1).state.min(:,2) = during_state_min;
limits(1).state.min(:,3) = final_state_min;

limits(1).state.max(:,1) = initial_state;
limits(1).state.max(:,2) = during_state_max;
limits(1).state.max(:,3) = final_state_max;

limits(1).control.min = [Tmin urmin utmin uhmin]';
limits(1).control.max = [Tmax urmax utmax uhmax]';

limits(1).parameter.min = [];
limits(1).parameter.max = [];

limits(1).path.min = 0;
limits(1).path.max = 0;



limits(1).event.min = r_event'; %trying to put a 1% error tolerance
limits(1).event.max = r_event';
%limits(1).event.min = [0 0 0]';
%limits(1).event.max = [0 0 0]';
setup.limits = limits;

%% GUESS YAGO

%initialguess = load("2.Initial_guesses/Moscow_Arc_8/SF_MBH/Results_Method_1/Moscow_Arc8_n20.mat").guess;
initialguess = load("2.Initial_guesses/Moscow_Full_Combined_MBH_MBH_Method_1").sol_guess;
guess_state_OE = zeros(length(initialguess.state(:,:,arc_index)),6);
guess_state_MEE = zeros(length(initialguess.state(:,:,arc_index)),7);
for j = [1:length(initialguess.state(:,:,arc_index))]
    [a,e,i,RAAN,w,MA,nu] = orbitalElements(initialguess.state(j,1:3,arc_index), initialguess.state(j,4:6,arc_index), pars);
    guess_state_OE(j,:) = [a,e,i,RAAN,w,nu];
    guess_state_MEE(j,1:6) = orbital2equinoctial(guess_state_OE(j,:));
    guess_state_MEE(j,7) = initialguess.mass(j,:,arc_index)./setup.m_init;
end
% guess_state_MEE(:,2) = guess_state_MEE(:,2)./fscale;
% guess_state_MEE(:,3) = guess_state_MEE(:,3)./gscale;
% guess_state_MEE(:,4) = guess_state_MEE(:,4)./hscale;
% guess_state_MEE(:,5) = guess_state_MEE(:,5)./kscale;
guess_state_MEE(:,6) = guess_state_MEE(:,6)./(4*pi);

r_vec = initialguess.state(:,1:3,arc_index);
v_vec = initialguess.state(:,4:6,arc_index);

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

t1 = initialguess.thrust(:,1,arc_index);
t2 = initialguess.thrust(:,2,arc_index);
t3 = initialguess.thrust(:,3,arc_index);
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




guess_control = zeros(length(initialguess.state(:,:,arc_index)),4);
guess_control(:,1) = magT/1000/(setup.m_init*setup.ac);
for i = 1:length(uMEE(1,1,:))
guess_control(i,2:4) = uMEE(:,:,i)';
end

%initialguess.parameter(1:length(initialguess.phase.control(:,1)),1) = initialguess.parameter;

guess(1).time = initialguess.times(:,:,arc_index)./setup.tc;
guess(1).state = guess_state_MEE;
guess(1).control = guess_control;
% guess(1).parameter = initialguess.parameter;

% plot3(initialguess.state(:,1), initialguess.state(:,2),initialguess.state(:,3))
% hold on
% quiver3(initialguess.state(:,1), initialguess.state(:,2), initialguess.state(:,3),initialguess.thrust(:,1),initialguess.thrust(:,2),initialguess.thrust(:,3))
% grid on
% plot3(arc_rf(1),arc_rf(2),arc_rf(3),'x')
% plot3(arc_ri(1),arc_ri(2),arc_ri(3),'o')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')

p = guess_state_MEE(:,1);
f = guess_state_MEE(:,2);%*fscale;
g = guess_state_MEE(:,3);%*gscale;
h = guess_state_MEE(:,4);%*hscale;
k = guess_state_MEE(:,5);%*kscale;
L = guess_state_MEE(:,6)*4*pi;

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

%% SECOND GUESS
prev_sol = load("MEE/guess_arc_direction_20_nodes_11000_iterations.mat").solution;
guess(1).time = prev_sol.time;
guess(1).state = prev_sol.state;
guess(1).control = prev_sol.control;
%%
setup.name = 'Single Arc Optimization MEE scaled';
setup.funcs.dae = 'MEEDaeScaled';
setup.funcs.cost = 'Cost';
setup.funcs.event = 'Event';
setup.auxdata = auxdata;
setup.guess = guess;
setup.linkages = [];
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


%%

time = solution.time;
p = solution.state(:,1);
f = solution.state(:,2);
g = solution.state(:,3);
h = solution.state(:,4);
k = solution.state(:,5);
L = solution.state(:,6);
w = solution.state(:,7);
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;
X = (r./s2).*(cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L));
Y = (r./s2).*(sin(L)-alpha2.*sin(L)+2.*h.*k.*cos(L));
Z = (2.*r./s2).*(h.*sin(L)-k.*cos(L));
vX = -(1./s2).*sqrt(1./p).*(sin(L)+alpha2.*sin(L)-2.*h.*k.*cos(L)+g-2.*f.*h.*k+alpha2.*g);
vY = -(1./s2).*sqrt(1./p).*(-cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L)-f+2.*g.*h.*k+alpha2.*f);
vZ = (2./s2).*sqrt(1./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);

r_vec = [X Y Z]; v_vec = [vX vY vZ];

hold on

plot3(X,Y,Z,'LineWidth',1.5)
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
grid on
view(3)
xf = solution.state(end,:);
pf = xf(1);
ff = xf(2);
gf = xf(3);
hf = xf(4);
kf = xf(5);
Lf = xf(6);
wf = xf(7);
q = 1+ff*cos(Lf)+gf*sin(Lf);
r = pf/q;
alpha2 = hf*hf-kf*kf;
chi = sqrt(hf*hf+kf*kf);
s2 = 1+chi*chi;

rX = (r/s2)*(cos(Lf)+alpha2*cos(Lf)+2*hf*kf*sin(Lf));
rY = (r/s2)*(sin(Lf)-alpha2*sin(Lf)+2*hf*kf*cos(Lf));
rZ = (2*r/s2)*(hf*sin(Lf)-kf*cos(Lf));
rVec = [rX rY rZ];

vX = -(1/s2)*sqrt(1/pf)*(sin(Lf)+alpha2*sin(Lf)-2*hf*kf*cos(Lf)+gf-2*ff*hf*kf+alpha2*gf);
vY = -(1/s2)*sqrt(1/pf)*(-cos(Lf)+alpha2*cos(Lf)+2*hf*kf*sin(Lf)-ff+2*gf*hf*kf+alpha2*ff);
vZ = (2/s2)*sqrt(1/pf)*(hf*cos(Lf)+kf*sin(Lf)+ff*hf+gf*kf);
vVec = [vX vY vZ]

% quiver3(rVec(1),rVec(2),rVec(3),vVec(1)*1e6,vVec(2)*1e6,vVec(3)*1e6)
VI_next = table2array(arcsTable(arc_index+1,4))/setup.vc;
deltaV = norm(vVec - VI_next)
% quiver3(rVec(1),rVec(2),rVec(3),VI_next(1)*1e6,VI_next(2)*1e6,VI_next(3)*1e6)


controls = solution.control;
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


for i = [1:height(solution.control)]
    uMEE(:,:,i) = solution.control(i,2:4);
end
%reverse to double check
for i = [1:length(Q_ECI2MEE(1,1,:))]
    uECI(:,:,i) = Q_ECI2MEE(:,:,i)'*uMEE(:,:,i);
end
uECIplot = zeros(height(solution.control),3);
for i = [1:height(solution.control)]

uECIplot(i,:) = uECI(:,:,i);
if solution.control(i,1) <= 0.0001;
    uECIplot(i,:) = [0 0 0];
end
end

quiver3(X,Y,Z,uECIplot(:,1),uECIplot(:,2),uECIplot(:,3))

%%
figure()

    hold on; grid on
    plot(solution.time*setup.tc/86400, solution.control(:,1)*1000*(setup.m_init*setup.ac),'LineWidth',1.5);
    xline((arc_tof)*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex');
    

xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Thrust$ (N)', 'FontSize', 14, 'Interpreter', 'latex');
title('Thrust History', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';
%% PLOT MASS HISTORY
figure()

    hold on; grid on
    plot(solution.time*setup.tc/86400, solution.state(:,end)*1500,'LineWidth',1.5);
    xline((arc_tof)*setup.tc/86400, '--', 'Color', 'r', 'LineWidth', 1, 'Interpreter', 'latex');
    

xlabel('$Time$ ([days])', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$Mass$ (kg)', 'FontSize', 14, 'Interpreter', 'latex');
title('Mass History', 'FontSize', 16, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.GridLineStyle = '--';