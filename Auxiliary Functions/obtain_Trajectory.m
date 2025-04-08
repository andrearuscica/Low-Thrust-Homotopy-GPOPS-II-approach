function trajectory = obtain_Trajectory(r_ini,v_ini,delta_t,pars,mesh_refiner)

%{
This function propagates the orbit providing the position vector of a
number of intermediate points in between the initial and final position:

r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
delta_t ---> Time of propagation [s]
pars ---> Structure containing problem data including mu
mesh_refiner ---> Number of points where the position is obtained

trajectory ---> Position vectors of the points [mesh_refiner x 3]
%}

%Obtain orbital elements of the segment
[a, e, i, RAN, Omega, MA_i] = Build_OE(r_ini, v_ini, pars.mu_sun);

%Obtain final position and velocity of the segment
[r_fin,v_fin] = propagate_kepler(r_ini,v_ini,delta_t,pars);
 
%Obtain initial and final true anomalies
TA_i = position2anomaly(r_ini, v_ini, pars.mu_sun, 0);
TA_f = position2anomaly(r_fin, v_fin, pars.mu_sun, 0);
if TA_i>TA_f
    TA_i = TA_i-360;
end

%Discretize true anomaly
TA_x = linspace(deg2rad(TA_i),deg2rad(TA_f),mesh_refiner);

%Apply the conic equation
r_x = zeros(size(TA_x));
for j = 1:length(TA_x)
    r_x(j) = a*(1-e^2)/(1+e*cos(TA_x(j)));
end

%Obtain the trajectory    
inc  = deg2rad(i);    %[rad] Body Inclination
RAAN = deg2rad(RAN);    %[rad] Body Longitude of Ascending Node
w    = deg2rad(Omega);    %[rad] Body Argument of Periapsis

trajectory = zeros(mesh_refiner,3);
for j = 1:length(r_x)
    rx = r_x(j)*(cos(TA_x(j) + w)*cos(RAAN) - sin(TA_x(j) + w)*cos(inc)*sin(RAAN));    %[km]
    ry = r_x(j)*(cos(TA_x(j) + w)*sin(RAAN) + sin(TA_x(j) + w)*cos(inc)*cos(RAAN));    %[km]
    rz = r_x(j)*(sin(TA_x(j) + w)*sin(inc));                                     %[km]
    
    trajectory(j,:) = [rx, ry, rz];
end