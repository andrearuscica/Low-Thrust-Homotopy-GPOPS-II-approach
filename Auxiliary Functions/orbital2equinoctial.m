function equinoctial = orbital2equinoctial(orbitalEl)
%   This functions computed the modified equinoctial elements from the
%   orbital elements.
%
%   INPUTS:
%       - orbitalEL: 1x6 vector containing the classical orbital elements [a [AU],e [-],i [°],RAAN [°],omega [°],nu [°]]
%   OUTPUTS:
%       - equinoctial: 1x6 vector containing the modified equinoctial
%       elements [p [AU], f [-], g [-], h [-], k [-], L [rad]]

a = orbitalEl(1);
e = orbitalEl(2);
i = deg2rad(orbitalEl(3));
RAAN = deg2rad(orbitalEl(4));
omega = deg2rad(orbitalEl(5));
nu = deg2rad(orbitalEl(6));

p = a*(1-e^2);
f = e*cos(omega + RAAN);
g = e*sin(omega + RAAN);
h = tan(i/2)*cos(RAAN);
k = tan(i/2)*sin(RAAN);
L = (RAAN + omega + nu);

equinoctial = [p, f, g, h, k, L];
end

