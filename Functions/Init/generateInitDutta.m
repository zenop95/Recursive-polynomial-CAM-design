function pp = generateInitDutta(orbit)

mu     = 3.986004418e14;    % [m^3/s^2]
HBR    = 4.295;                % [m]
a      = 6778000;        %[m] semimajor axis
e      = 0;              %[-] eccentricity
theta0 = pi;              %[rad] initial true anomaly
omega  = 0;              %[rad] argument of periapsis
RAAN   = 0;              %[rad] right ascension node longitude
inc    = 0.87266;              %[rad] inclination
n      = (mu/a^3)^(1/2);
T      = 2*pi/n;         %[s] orbital period

aD = 6753000;
eD = 0.0037;
theta0D = pi;
incD = 0;
omegaD = 0;
RAAND = 0;

%% Inertial state of the two objects
%inertial state of the satellite
x_s = kepler2cartesian(a,e,RAAN,inc,omega,theta0);
x_d = kepler2cartesian(aD,eD,RAAND,incD,omegaD,theta0D);

%% Relative initial conditions in ECI reference frame
C0  = [500^2*2*eye(3) zeros(3); zeros(3,6)];

pp = struct( ...
            'orbit', orbit, ...
            'mu',    mu, ...
            'HBR',   HBR, ...
            'x_s',   x_s, ...
            'x_d',   x_d, ...
            'C0',    C0, ...
            'C00',    C0, ...
            'T',     T, ...
            'n',     n ...
            );
end