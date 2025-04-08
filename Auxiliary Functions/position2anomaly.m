function Anomaly = position2anomaly(r_vec, v_vec, mu, select)

%{
This function allows the computation of the anomalies of an orbit. The
algorithm is prepared for elliptic and hyperbolic orbits. The inputs are:
[r_vec] ---> The position vector [L]
[v_vec] ---> The velocity vector [L/T]
[mu] ---> The gravitational parameter [L^3/T^2]
[select] ---> Several values are possible:
              0 - True anomaly [deg] (elliptic orbit: [0,360])
              1 - Eccentric anomaly [deg] (elliptic orbit: [0,360])
              2 - Mean anomaly [deg] (elliptic orbit: [0,360])
%}

%Compute angular momentum
h_vec = cross(r_vec,v_vec);
h = norm(h_vec);
if h<10e-2
    disp("Special case found: h=0")
end

%Compute eccentricity
e_vec = (cross(v_vec,h_vec)/mu) - (r_vec/norm(r_vec));
e = norm(e_vec);
if e<10e-2
    disp("Special case found: circular orbit")
elseif e ==1
    disp("Special case found: parabolic orbit")
end

vr = dot(v_vec,r_vec/norm(r_vec));
if e<1
    %Compute the true anomaly
    if vr >= 0 
        nu = acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    else
        nu = 2*pi - acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    end
    
    %Compute eccentric anomaly
    E = atan2(sqrt(1-e^2)*sin(nu),e+cos(nu));
    if E<0
        E = E + 2*pi;
    end
    
    %Compute the mean anomaly
    MA = rad2deg(E - e*sin(E));

elseif e>1
    %Compute the true anomaly
    if vr >= 0 
        nu = acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    else
        nu = -acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    end

    %Compute eccentric anomaly
    E = 2*atanh(sqrt((e-1)/(e+1))*tan(nu/2));

    %Compute the mean anomaly
    MA = rad2deg(e*sinh(E)-E);
end

if select == 0
    Anomaly = rad2deg(nu);
elseif select == 1
    Anomaly = rad2deg(E);
else
    Anomaly = MA;
end

