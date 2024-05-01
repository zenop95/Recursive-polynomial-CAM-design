function pp = generateInitGeo(orbit)

mu  = 398600.4418;    % [km^3/s^2]

%% Inertial state of the Primary
r0 = [9334.75561447905; 41120.9856893676; 1.57083557118839]; % From corrected Laura's case 1
v0 = [-2.99818420674816; .680764971686596; .00637727006252368];
cart = [r0; v0];
p = cartesian2kepler(cart,mu);
p.cart0 = cart;
p.x0    = cart;
p.C0     = 0*diag(([4.5 10 5 0.15 0.75 0.25]/1e3).^2);  % [km^2] [km^2/s^2] Covariance at TCA
p.HBR     = 0.05;  % [km^2] [km^2/s^2] Covariance at TCA
p.mass   = 500;            % [kg] mass
p.A_drag = 1;              % [m^2] drag surface area
p.Cd     = 2.2;            % [-] shape coefficient for drag
p.A_srp  = 1;              % [m^2] SRP surface area
p.Cr     = 1.31;           % [-] shape coefficient for SRP
T = 2*pi/p.n;
%% Relative initial conditions in Hill reference frame
% A0      = 0.500;
% B0      = 0.100;
% y_off   = 0;
% c       = 5;
% alpha0  = pi/180*15*(c-3);
% beta0   = pi/2 + alpha0;
% n       = p.n;
% dr0      = [A0*cos(alpha0); -2*A0*sin(alpha0) + y_off; B0*cos(beta0)];     % [km] Relative position at TCA in LVLH
% dv0      = [-n*A0*sin(alpha0); -2*n*A0*cos(alpha0); -n*B0*sin(beta0)];     % [km/s] Relative velocity at TCA in LVLH
dr0      = [0; 0.5; 0];     % [km] Relative position at TCA in LVLH
dv0      = [-1e-3; 0; 0];     % [km/s] Relative velocity at TCA in LVLH

%% Relative state of the secondary object at TCA
j = 1;
s(j).tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
s(j).x0       = [];
s(j).relState = [dr0; dv0];                                                % [km] [km/s] Relative cartesian state at TCA in RTN
s(j).C0       = diag([.45 1 .05 0.0015 0.0075 0.0025]).^2;                    % [km^2] [km^2/s^2] Covariance at TCA
s(j).HBR      = p.HBR + 0.033;  % [km]
s(j).mass     = 500;          % [kg] mass
s(j).A_drag   = 1;            % [m^2] drag surface area
s(j).Cd       = 2.2;          % [-] shape coefficient for drag
s(j).A_srp    = 1;            % [m^2] SRP surface area
s(j).Cr       = 1.31;         % [-] shape coefficient for SRP
s(j).x          = [];         % [-] shape coefficient for SRP
s(j).covariance = [];         % [-] shape coefficient for SRP
s(j).w     = 1;
s(j).cdm   = true;
pp = struct( ...
            'orbit',     orbit, ...
            'mu',        mu, ...
            'primary',   p, ...
            'secondary', s, ...
            'T',         T ...
            );
end