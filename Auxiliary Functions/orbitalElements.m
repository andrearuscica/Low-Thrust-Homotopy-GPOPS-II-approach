function [a,e,i,RAAN,w,MA,nu] = orbitalElements(r_vec, v_vec, pars)


% Compute the main orbital elements from radius and velocity vectors
% INPUTS:
% - r_vec: radius vector [km];
% - v_vec: velocity vector [km/s];
% OUTPUT:
% [a,e,i,RAAN,w,MA]: vector containing the orbital elements:
%       a: Semimajor axis [AU];
%       e: Eccentricity [-];
%       i: Inclination [°];
%       RAAN: Right ascension of the ascending node [°];
%       w: Argument of periapsis [°];
%       MA: Mean anomaly [°];


r = norm(r_vec);
v = norm(v_vec);
v_r = dot(r_vec/r, v_vec); % tangential velocity versor
%v_p = sqrt(v^2 - v_r^2);   %perpendicular velocity versor
 

% Orbital Angular Momentum
h_vec = cross(r_vec, v_vec);
h = norm(h_vec);

% Inclination
i = acos(h_vec(3) / h);
i = i*180/pi;

% Right ascension of the ascending node
K = [0, 0, 1];
N_vec = cross(K, h_vec);
N = norm(N_vec);
if N_vec(2) <= 0
    RAAN = 2 * pi - acos(N_vec(1)/N);
else
    RAAN = acos(N_vec(1)/N);
end

% Eccentricity
e_vec = cross(v_vec, h_vec) / pars.mu_sun - r_vec / r;
e = norm(e_vec);

% Argument of periapsis
w = acos(dot(N_vec, e_vec) / (N * e));
if e_vec(3) < 0
    w = 2 * pi - w; 
end


% True anomaly
nu = acos(dot(e_vec, r_vec) / (e * r));
if dot(r_vec, v_vec) < 0
    nu = 2 * pi - nu;
end


% Semi major axis
a = r * (1 + e * cos(nu)) / (1 - e^2);

% Mean anomaly
tanE_half = tan(nu/2) * sqrt((1-e)/(1+e));
E = 2* atan(tanE_half);
MA = E - e*sin(E);
if MA < 0
    MA = 2*pi + MA;
end

RAAN = RAAN*180/pi;
MA = MA*180/pi;
nu = nu*180/pi;
a = a/pars.AU;
w = w*180/pi;
end