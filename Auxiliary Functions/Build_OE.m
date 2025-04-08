function [a,e,i,RAAN,w,MA] = Build_OE(r_vec, v_vec,mu)

%{
This function builds the orbital elements of an object from its current
position and speed. It can't compute pure parabolic orbits.

INPUTS:
r_vec ---> Position vector in cartesian coordinates [L]
v_vec ---> Velocity vector in cartesian coordinates [L/T]
mu ---> Gravitational Parameter [L^3/T^2]

OUTPUTS:
a ---> Semi-major axis [L]
e ---> Eccentricity []
i ---> Inclination [deg] (0<i<180)
RAAN ---> Right ascension of the ascending node [deg] (0<RAAN<360)
w ---> Longitude of the periapsis [deg] (0<w<360)
MA ---> Mean anomaly [deg] (elliptic: 0<MA<360, hyperbolic: -inf<MA<inf)

%}

%Compute angular momentum
h_vec = cross(r_vec,v_vec);
h = norm(h_vec);
if h<10e-4
    warning("Special case found: h=0")
end

%Compute eccentricity
e_vec = (cross(v_vec,h_vec)/mu) - (r_vec/norm(r_vec));
e = norm(e_vec);
if e<10e-3
    warning("Special case found: circular orbit")
elseif e ==1
    warning("Special case found: parabolic orbit")
end

%Compute energy
Energy = ((norm(v_vec)^2)/2)-(mu/norm(r_vec));

%Compute semi-major axis
a = -mu/(2*Energy);

%Compute inclination
hz = h_vec(3);
i = rad2deg(acos(hz/h));
if i < 10e-2
    warning("Special case found: i=0")
end

%Compute the RAAN
n = cross([0 0 1],h_vec);
if n(2)>= 0
    RAAN = rad2deg(acos(n(1)/norm(n)));
else
    RAAN = rad2deg(2*pi - acos(n(1)/norm(n)));
end

%Compute the argument of the periapsis
if e_vec(3)>=0
    w = rad2deg(acos(dot(n,e_vec)/(norm(n)*norm(e_vec))));
else
    w = rad2deg(2*pi - acos(dot(n,e_vec)/(norm(n)*norm(e_vec))));
end

%Compute the true anomaly
vr = dot(v_vec,r_vec/norm(r_vec));
if e<1
    if vr >= 0 
        nu = acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    else
        nu = 2*pi - acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    end
elseif e>1
    if vr >= 0 
        nu = acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    else
        nu = -acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    end
end

%Compute eccentric anomaly
if e<1
    E = atan2(sqrt(1-e^2)*sin(nu),e+cos(nu));
    if E<0
        E = E + 2*pi;
    end
elseif e>1
    H = 2*atanh(sqrt((e-1)/(e+1))*tan(nu/2));
end

%Compute the mean anomaly
if e<1
    MA = rad2deg(E - e*sin(E));
elseif e>1
    MA = rad2deg(e*sinh(H)-H);
end


