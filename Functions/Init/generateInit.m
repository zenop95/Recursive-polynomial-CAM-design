function pp = generateInit(orbit)

mu  = 398600.4418;    % [km^3/s^2]

%% Inertial State of the Primary spacecraft
p.C0      = .1*diag([.005 .01 .0045 2.5e-4 7.5e-4 1.5e-4]).^2;             % [km^2] [km^2/s^2] Covariance at TCA
p.ecc    = 0;              % [-] eccentricity 
p.theta0 = 0;              % [rad] initial true anomaly
p.omega  = 0;              % [rad] argument of periapsis
p.RAAN   = 0;              % [rad] right ascension node longitude
p.inc    = deg2rad(53);    % [rad] inclination
p.a      = 6928;           % [km] semimajor axis
p.n      = (mu/p.a^3)^(1/2);    %[rad/s] mean motion
T        = 2*pi/p.n;         % [s] orbital period
p.HBR    = 0.003;           % [km]
p.mass   = 260;            % [kg] mass
p.A_drag = 1;              % [m^2] drag surface area
p.Cd     = 2.2;            % [-] shape coefficient for drag
p.A_srp  = 1;              % [m^2] SRP surface area
p.Cr     = 1.31;           % [-] shape coefficient for SRP
p.cart0  = kepler2cartesian(p.a,p.ecc,p.RAAN,p.inc,p.omega,p.theta0,mu);
% p.inc    = 1.6;            % [rad] inclination
% p.n      = 1.125915e-3;    %[rad/s] mean motion
% p.a      = (mu/p.n^2)^(1/3); % [km] semimajor axis
% T        = 2*pi/p.n;         % [s] orbital period
% p.HBR    = 0.02;           % [km]
% p.mass   = 500;            % [kg] mass
% p.A_drag = 1;              % [m^2] drag surface area
% p.Cd     = 2.2;            % [-] shape coefficient for drag
% p.A_srp  = 1;              % [m^2] SRP surface area
% p.Cr     = 1.31;           % [-] shape coefficient for SRP

%% Relative initial conditions in Hill reference frame
A0      = 1.000;
B0      = 0.010;
y_off   = 0.1;
c       = 5;
alpha0  = pi/180*15*(c-3);
beta0   = pi/2 + alpha0;
n       = p.n;
r0      = [A0*cos(alpha0); -2*A0*sin(alpha0) + y_off; B0*cos(beta0)];      % [km] Relative position at TCA in LVLH
v0      = [-n*A0*sin(alpha0); -2*n*A0*cos(alpha0); -n*B0*sin(beta0)];      % [km/s] Relative velocity at TCA in LVLH

%% Relative state of the secondary object at TCA
j = 1;
s(j).tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
s(j).x0       = [];
% s(j).relState = -[r0; v0];                                                  % [km] [km/s] Relative cartesian state at TCA in RTN
s(j).relState = [1 0.1 0 0 -2e-3 -1e-5]';                        % [km] [km/s] Relative cartesian state at TCA
s(j).C0       = diag([.05 .1 .045 2.5e-4 7.5e-4 1.5e-4]).^2;             % [km^2] [km^2/s^2] Covariance at TCA
s(j).HBR      = p.HBR + 0.032;  % [km]
s(j).mass     = 200;          % [kg] mass
s(j).A_drag   = 1;            % [m^2] drag surface area
s(j).Cd       = 2.2;          % [-] shape coefficient for drag
s(j).A_srp    = 1;            % [m^2] SRP surface area
s(j).Cr       = 1.31;         % [-] shape coefficient for SRP
s(j).x          = [];         % [-] shape coefficient for SRP
s(j).covariance = [];         % [-] shape coefficient for SRP
s(j).w          = 1;
s(j).cdm        = true;
s(j).ang        = false;
pp = struct( ...
            'orbit',     orbit, ...
            'mu',        mu, ...
            'primary',   p, ...
            'secondary', s, ...
            'T',         T ...
            );
end