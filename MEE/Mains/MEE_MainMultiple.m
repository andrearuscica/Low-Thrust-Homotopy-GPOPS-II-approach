clear all; clc; %close all;


nphases = 2;
starting_arc = 3;

arc_index = 12;
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

%arcs
arcsTable = problem_data.arcs; %i for rows, ri = 2,rf=3,vi = 4, vf= 5, tof = 6
auxdata.arcsTable = arcsTable;
arc_ri = table2array(arcsTable(arc_index,2));
arc_rf = table2array(arcsTable(arc_index,3));
arc_vi = table2array(arcsTable(arc_index,4));
arc_tof = table2array(arcsTable(arc_index,6));
arc_prev_vf = table2array(arcsTable(arc_index-1,5));
w0 = table2array(arcsTable(starting_arc,7));

auxdata.w0 = w0;
t0 = 0; tmin = t0; tf = arc_tof; tfmin = 0.6*tf; tmax = 1.25*tf;

MEE_Table = problem_data.MEE;

previous_arc_MEE = table2array(MEE_Table(starting_arc-1,2));
p0 = previous_arc_MEE(1)*AU;
f0 = previous_arc_MEE(2);
g0 = previous_arc_MEE(3);
h0 = previous_arc_MEE(4);
k0 = previous_arc_MEE(5);
L0 = previous_arc_MEE(6);

fmin = -1; fmax = +1;
gmin = -1; gmax = +1;
hmin = -1; hmax = +1;
kmin = -1; kmax = +1;
Lmin = L0; Lmax = 9*2*pi;
wmin = 500; wmax = w0;
urmin = -1; urmax = +1;
utmin = -1; utmax = +1;
uhmin = -1; uhmax = +1;
Tmin = 0; Tmax = T;

initial_state = [p0 f0 g0 h0 k0 L0 w0]';
during_state_min = [0.3*p0 fmin gmin hmin kmin Lmin wmin]';
during_state_max = [3*p0 fmax gmax hmax kmax Lmax wmax]'; 
final_state_min = [0.3*p0 fmin gmin hmin kmin Lmin wmin]';
final_state_max = [3*p0 fmax gmax hmax kmax Lmax wmax]';
% initial MEE dell'orbita del previous lambert;
% final a piacere;



%% LIMITS

for iphase = [1:nphases]
    iphase
    num_nodes = 10; limits(iphase).nodes = num_nodes;

    arc_tof = table2array(arcsTable(iphase+starting_arc-1,6));

    tf = arc_tof; tfmin = 0.6*tf; tmax = 2*tf;
    
    arc_ri = table2array(arcsTable(iphase+starting_arc-1,2)); % +2 because I start from arc 3
    arc_rf = table2array(arcsTable(iphase+starting_arc-1,3));
    arc_vi = table2array(arcsTable(iphase+starting_arc-1,4));
    arc_tof = table2array(arcsTable(iphase+starting_arc-1,6));
    arc_prev_vf = table2array(arcsTable(iphase+starting_arc-2,5));
    arc_VI_next = table2array(arcsTable(iphase+starting_arc,4));
    

    
    

limits(iphase).time.min = [0 tfmin];
limits(iphase).time.max = [0 tmax];


    

    if iphase == 1
    limits(iphase).state.min(:,1) = initial_state;
    limits(iphase).state.min(:,2) = during_state_min;
    limits(iphase).state.min(:,3) = final_state_min;
    
    limits(iphase).state.max(:,1) = initial_state;
    limits(iphase).state.max(:,2) = during_state_max;
    limits(iphase).state.max(:,3) = final_state_max;
    else
    limits(iphase).state.min(:,1) = during_state_min;
    limits(iphase).state.min(:,2) = during_state_min;
    limits(iphase).state.min(:,3) = final_state_min;
        
    limits(iphase).state.max(:,1) = during_state_max;
    limits(iphase).state.max(:,2) = during_state_max;
    limits(iphase).state.max(:,3) = final_state_max;

    end
    
    limits(iphase).control.min = [Tmin urmin utmin uhmin]';
    limits(iphase).control.max = [Tmax*10 urmax utmax uhmax]';
    
    limits(iphase).parameter.min = [];
    limits(iphase).parameter.max = [];
    
    limits(iphase).path.min = [0];
    limits(iphase).path.max = [0];
   
    % limits(iphase).event.min = [0 0 0]';
    % limits(iphase).event.max = [0 0 0]';
    limits(iphase).event.min = arc_rf';
    limits(iphase).event.max = arc_rf';
    
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
% setup.limits = limits;
% num_nodes = 20; limits(1).nodes = num_nodes;
% 
% limits(1).time.min = [0 tfmin];
% limits(1).time.max = [0 tmax];
% 
% limits(1).state.min(:,1) = initial_state;
% limits(1).state.min(:,2) = during_state_min;
% limits(1).state.min(:,3) = final_state_min;
% 
% limits(1).state.max(:,1) = initial_state;
% limits(1).state.max(:,2) = during_state_max;
% limits(1).state.max(:,3) = final_state_max;
% 
% limits(1).control.min = [Tmin urmin utmin uhmin]';
% limits(1).control.max = [Tmax urmax utmax uhmax]';
% 
% limits(1).parameter.min = [];
% limits(1).parameter.max = [];
% 
% limits(1).path.min = 0;
% limits(1).path.max = 0;
% 
% %limits(1).event(min-max) = [rasteroide]
% limits(1).event.min = r_event';
% limits(1).event.max = r_event';
setup.limits = limits;

%% GUESS Constant thrust
%SCALING
setup.lc  = 149597870.700;
setup.mu  = 132712440018;
setup.tc  = sqrt(setup.lc^3/setup.mu);
setup.vc  = setup.lc/setup.tc;
setup.ac  = setup.lc/setup.tc^2;
setup.m_init = 1500;

guessData.Isp = 3000;
guessData.T = T*10;

ODEOPTS = odeset('RelTol',1e-6,'AbsTol',1e-6);
guesses = zeros(nphases,7);
for i = [1:nphases]
    i
    previous_arc_MEE = table2array(MEE_Table(i+starting_arc-2,2));
    pg = previous_arc_MEE(1);%*AU;
    fg = previous_arc_MEE(2);
    gg = previous_arc_MEE(3);
    hg = previous_arc_MEE(4);
    kg = previous_arc_MEE(5);
    Lg = previous_arc_MEE(6);
    arc_ri = table2array(arcsTable(i+starting_arc-1,2));
    arc_prev_vf = table2array(arcsTable(i+starting_arc-2,5));
    arc_rf = table2array(arcsTable(i+starting_arc-1,3));
    arc_tof = table2array(arcsTable(i+starting_arc-1,6));
    wg = table2array(arcsTable(i+starting_arc-1,7));

    dir1 = arc_rf - arc_ri;
    guessData.direction = dir1./norm(dir1);
    tic;
    xx = [pg fg gg hg kg Lg wg]';
    [TOUT, XOUT] = ode45(@MEE_FirstGuessScaledProva, [0,arc_tof./setup.tc],xx,ODEOPTS,guessData);
    time = toc
    guess(i).state = XOUT;
    guess(i).time = TOUT;

    p_guess = XOUT(:,1);
    f_guess = XOUT(:,2);
    g_guess = XOUT(:,3);
    h_guess = XOUT(:,4);
    k_guess = XOUT(:,5);
    L_guess = XOUT(:,6);
    w_guess = XOUT(:,7);


q = 1+f_guess.*cos(L_guess)+g_guess.*sin(L_guess);
r = p_guess./q;
alpha2 = h_guess.*h_guess-k_guess.*k_guess;
chi = sqrt(h_guess.*h_guess+k_guess.*k_guess);
s2 = 1+chi.*chi;

% ----------------------------------------------------------------------- %
% Trajectory in Cartesian
% ----------------------------------------------------------------------- % 
X = (r./s2).*(cos(L_guess)+alpha2.*cos(L_guess)+2.*h_guess.*k_guess.*sin(L_guess));
Y = (r./s2).*(sin(L_guess)-alpha2.*sin(L_guess)+2.*h_guess.*k_guess.*cos(L_guess));
Z = (2.*r./s2).*(h_guess.*sin(L_guess)-k_guess.*cos(L_guess));
vX = -(1./s2).*sqrt(mu./p_guess).*(sin(L_guess)+alpha2.*sin(L_guess)-2.*h_guess.*k_guess.*cos(L_guess)+g_guess-2.*f_guess.*h_guess.*k_guess+alpha2.*g_guess);
vY = -(1./s2).*sqrt(mu./p_guess).*(-cos(L_guess)+alpha2.*cos(L_guess)+2.*h_guess.*k_guess.*sin(L_guess)-f_guess+2.*g_guess.*h_guess.*k_guess+alpha2.*f_guess);
vZ = (2./s2).*sqrt(mu./p_guess).*(h_guess.*cos(L_guess)+k_guess.*sin(L_guess)+f_guess.*h_guess+g_guess.*k_guess);

plot3(X,Y,Z)
grid on; hold on;
%plot3(arc_rf(1),arc_rf(2),arc_rf(3),'x')
%plot3(arc_ri(1),arc_ri(2),arc_ri(3),'o')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')

r_vec = [X Y Z];
v_vec = [vX vY vZ];

normR = (r_vec(:,1).^2 + r_vec(:,2).^2 + r_vec(:,3).^2).^(1/2);
normV = sqrt(v_vec(:,1).^2 + v_vec(:,2).^2 + v_vec(:,3).^2).^(1/2);

ir = r_vec./normR;
h_vec = cross(r_vec,v_vec,2);
normH = (h_vec(:,1).^2 + h_vec(:,2).^2 + h_vec(:,3).^2).^(1/2);
in = (h_vec)./normH;
it = cross(in,ir,2);

    for j = [1:height(ir)]
        %each "page" of the matrix is a conversion matrix for a state
        Q_ECI2MEE(:,:,j) = [ir(j,:); it(j,:); in(j,:)];
    end

    for j = [1:length(Q_ECI2MEE(1,1,:))]
        uMEE(:,:,j) = Q_ECI2MEE(:,:,j)*guessData.direction';
    end

guess_control = zeros(length(XOUT),4);
guess_control(:,1) = T;
    for j = 1:length(uMEE(1,1,:))
    guess_control(j,2:4) = uMEE(:,:,j)';
    end


guess(i).control = guess_control;
guess(i).parameter = [];


hold on
%plot3(XOUT(:,1),XOUT(:,2),XOUT(:,3))

%plot3(arc_rf(1),arc_rf(2),arc_rf(3),'x')
%plot3(arc_ri(1),arc_ri(2),arc_ri(3),'o')
%quiver3(arc_ri(1),arc_ri(2),arc_ri(3),dir1(1),dir1(2),dir1(3))
view(3)

end






%% GUESS YAGO

initialguess = load("2.Initial_guesses/Moscow_Arc_8/SF_MBH/Results_Method_1/Moscow_Arc8_n20.mat").guess;
guess_state_OE = zeros(length(initialguess.state),6);
guess_state_MEE = zeros(length(initialguess.state),7);
for j = [1:length(initialguess.state)]
    [a,e,i,RAAN,w,MA,nu] = orbitalElements(initialguess.state(j,1:3), initialguess.state(j,4:6), pars);
    guess_state_OE(j,:) = [a,e,i,RAAN,w,nu];
    guess_state_MEE(j,1:6) = orbital2equinoctial(guess_state_OE(j,:));
    guess_state_MEE(j,7) = initialguess.mass(j);
end
guess_state_MEE(:,1) = guess_state_MEE(:,1)*AU;

r_vec = initialguess.state(:,1:3);
v_vec = initialguess.state(:,4:6);

normR = (r_vec(:,1).^2 + r_vec(:,2).^2 + r_vec(:,3).^2).^(1/2);
normV = sqrt(v_vec(:,1).^2 + v_vec(:,2).^2 + v_vec(:,3).^2).^(1/2);

ir = r_vec./normR;
h_vec = cross(r_vec,v_vec,2);
normH = (h_vec(:,1).^2 + h_vec(:,2).^2 + h_vec(:,3).^2).^(1/2);
in = (h_vec)./normH;
it = cross(in,ir,2);

for i = [1:height(ir)]
    %each "page" of the matrix is a conversion matrix for a state
    Q_ECI2MEE(:,:,i) = [ir(i,:); it(i,:); in(i,:)];
end

t1 = initialguess.thrust(:,1);
t2 = initialguess.thrust(:,2);
t3 = initialguess.thrust(:,3);
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




guess_control = zeros(length(initialguess.state),4);
guess_control(:,1) = magT;
for i = 1:length(uMEE(1,1,:))
guess_control(i,2:4) = uMEE(:,:,i)';
end

%initialguess.parameter(1:length(initialguess.phase.control(:,1)),1) = initialguess.parameter;

guess(1).time = initialguess.times;
guess(1).state = guess_state_MEE;
guess(1).control = guess_control;
% guess(1).parameter = initialguess.parameter;

plot3(initialguess.state(:,1), initialguess.state(:,2),initialguess.state(:,3))
hold on
quiver3(initialguess.state(:,1), initialguess.state(:,2), initialguess.state(:,3),initialguess.thrust(:,1),initialguess.thrust(:,2),initialguess.thrust(:,3))
grid on
plot3(arc_rf(1),arc_rf(2),arc_rf(3),'x')
plot3(arc_ri(1),arc_ri(2),arc_ri(3),'o')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')

p = guess_state_MEE(:,1);
f = guess_state_MEE(:,2);
g = guess_state_MEE(:,3);
h = guess_state_MEE(:,4);
k = guess_state_MEE(:,5);
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

%plot3(X,Y,Z)


%%
setup.name = 'Multiple Arc Optimization MEE';
setup.funcs.dae = 'Dae';
setup.funcs.cost = 'MEE_CostMultiple';
setup.funcs.event = 'Event';
setup.auxdata = auxdata;
setup.guess = guess;
setup.linkages = linkages;
setup.funcs.link  = 'CartLinks';
setup.derivatives = 'automatic';
setup.parallel    = 'no';
setup.autoscale = 'on'; %Remember to on for pseudospectral method
setup.solver ='ipopt';
setup.method ='pseudospectral'; %hermite-simpson

 
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
p = solution(1).state(:,1);
f = solution(1).state(:,2);
g = solution(1).state(:,3);
h = solution(1).state(:,4);
k = solution(1).state(:,5);
L = solution(1).state(:,6);
w = solution(1).state(:,7);
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;
X = (r./s2).*(cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L));
Y = (r./s2).*(sin(L)-alpha2.*sin(L)+2.*h.*k.*cos(L));
Z = (2.*r./s2).*(h.*sin(L)-k.*cos(L));
vX = -(1./s2).*sqrt(mu./p).*(sin(L)+alpha2.*sin(L)-2.*h.*k.*cos(L)+g-2.*f.*h.*k+alpha2.*g);
vY = -(1./s2).*sqrt(mu./p).*(-cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L)-f+2.*g.*h.*k+alpha2.*f);
vZ = (2./s2).*sqrt(mu./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);

r_vec = [X Y Z]; v_vec = [vX vY vZ];
hold on
view(2)
plot3(X,Y,Z)


p = solution(2).state(:,1);
f = solution(2).state(:,2);
g = solution(2).state(:,3);
h = solution(2).state(:,4);
k = solution(2).state(:,5);
L = solution(2).state(:,6);
w = solution(2).state(:,7);
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;
X = (r./s2).*(cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L));
Y = (r./s2).*(sin(L)-alpha2.*sin(L)+2.*h.*k.*cos(L));
Z = (2.*r./s2).*(h.*sin(L)-k.*cos(L));
vX = -(1./s2).*sqrt(mu./p).*(sin(L)+alpha2.*sin(L)-2.*h.*k.*cos(L)+g-2.*f.*h.*k+alpha2.*g);
vY = -(1./s2).*sqrt(mu./p).*(-cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L)-f+2.*g.*h.*k+alpha2.*f);
vZ = (2./s2).*sqrt(mu./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);

r_vec = [X Y Z]; v_vec = [vX vY vZ];
hold on
view(2)
plot3(X,Y,Z)



%%
for phase = [1:2]
xf = solution(phase).state(end,:);
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

vX = -(1/s2)*sqrt(mu/pf)*(sin(Lf)+alpha2*sin(Lf)-2*hf*kf*cos(Lf)+gf-2*ff*hf*kf+alpha2*gf);
vY = -(1/s2)*sqrt(mu/pf)*(-cos(Lf)+alpha2*cos(Lf)+2*hf*kf*sin(Lf)-ff+2*gf*hf*kf+alpha2*ff);
vZ = (2/s2)*sqrt(mu/pf)*(hf*cos(Lf)+kf*sin(Lf)+ff*hf+gf*kf);
vVec = [vX vY vZ];

quiver3(rVec(1),rVec(2),rVec(3),vVec(1)*1e6,vVec(2)*1e6,vVec(3)*1e6)
VI_next = table2array(arcsTable(arc_index+1,4));
deltaV = norm(vVec - VI_next)
quiver3(rVec(1),rVec(2),rVec(3),VI_next(1)*1e6,VI_next(2)*1e6,VI_next(3)*1e6)


controls = solution(phase).control;
normR = (r_vec(:,1).^2 + r_vec(:,2).^2 + r_vec(:,3).^2).^(1/2);
normV = sqrt(v_vec(:,1).^2 + v_vec(:,2).^2 + v_vec(:,3).^2).^(1/2);

ir = r_vec./normR;
h_vec = cross(r_vec,v_vec,2);
normH = (h_vec(:,1).^2 + h_vec(:,2).^2 + h_vec(:,3).^2).^(1/2);
in = (h_vec)./normH;
it = cross(in,ir,2);

for i = [1:height(ir)]
    %each "page" of the matrix is a conversion matrix for a state
    Q_ECI2MEE(:,:,i) = [ir(i,:); it(i,:); in(i,:)];
end


for i = [1:height(solution(phase).control)]
    uMEE(:,:,i) = solution(phase).control(i,2:4);
end
%reverse to double check
for i = [1:length(Q_ECI2MEE(1,1,:))]
    uECI(:,:,i) = Q_ECI2MEE(:,:,i)'*uMEE(:,:,i);
end
uECIplot = zeros(height(solution(phase).control),3);
for i = [1:height(solution(phase).control)]
uECIplot(i,:) = uECI(:,:,i);
end
quiver3(X,Y,Z,uECIplot(:,1),uECIplot(:,2),uECIplot(:,3))

end