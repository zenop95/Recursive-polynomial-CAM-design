function pp = generateInitMasson(orbit)

mu  = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = diag(([4.5 10 5 0.15 0.75 0.25]/1e3).^2)/4;  % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = .00145;              % [-] eccentricity 
primary.theta0 = pi/2;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = deg2rad(86.4);              % [rad] inclination
primary.n      = 1e-3;           %[rad/s] mean motion
primary.a      = (mu/primary.n^2)^(1/3); % [km] semimajor axis
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = 0.02;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0 = kepler2cartesian(primary.a,primary.e,primary.RAAN,primary.inc,primary.omega,primary.theta0,mu);

%% Secondary structure
secondary = struct();
secondary(1).tca      = 2*T;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary(1).x0       = [];         % [-] shape coefficient for SRP
secondary(1).relState = [-.0079 .035 .08 .0057 -12.5526 5.489]';                          % [km] [km/s] Relative cartesian state at TCA
secondary(1).C0       = 3/4*diag(([4.5 10 5 0.15 0.75 0.25]/1e3).^2);          % [km^2] [km^2/s^2] Covariance at TCA
secondary(1).HBR      = primary.HBR + 0.01;         % [km]
secondary(1).mass     = 100;          % [kg] mass
secondary(1).A_drag   = 1;            % [m^2] drag surface area
secondary(1).Cd       = 2.2;          % [-] shape coefficient for drag
secondary(1).A_srp    = 1;            % [m^2] SRP surface area
secondary(1).Cr       = 1.31;         % [-] shape coefficient for SRP
secondary(1).x          = [];         % [-] shape coefficient for SRP
secondary(1).covariance = [];         % [-] shape coefficient for SRP

secondary(2).tca      = 2*T+1.5*3600;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary(2).x0          = [];         % [-] shape coefficient for SRP
secondary(2).relState = [-.0168 -.0018 .041 -.038 -14.9108 -.7855]';                          % [km] [km/s] Relative cartesian state at TCA
secondary(2).C0       = 3/4*diag(([4.5 10 5 0.15 0.75 0.25]/1e3).^2);          % [km^2] [km^2/s^2] argument of periapsis
secondary(2).HBR      = primary.HBR +  0.015;         % [km]
secondary(2).mass     = 100;          % [kg] mass
secondary(2).A_drag   = 1;            % [m^2] drag surface area
secondary(2).Cd       = 2.2;          % [-] shape coefficient for drag
secondary(2).A_srp    = 1;            % [m^2] SRP surface area
secondary(2).Cr       = 1.31;         % [-] shape coefficient for SRP
secondary(2).x          = [];         % [-] shape coefficient for SRP
secondary(2).covariance = [];         % [-] shape coefficient for SRP

secondary(3).tca      = 2*T+10*3600;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary(3).x0       = [];         % [-] shape coefficient for SRP
secondary(3).relState = [.0009 -.0352 -.0003 .0006 -.0005 .0058]';                          % [km] [km/s] Relative cartesian state at TCA
secondary(3).C0       = 3/4*diag(([4.5 10 5 0.15 0.75 0.35]/1e3).^2);          % [km^2] [km^2/s^2] argument of periapsis
secondary(3).HBR      = primary.HBR +  0.02;         % [km]
secondary(3).mass     = 200;          % [kg] mass
secondary(3).A_drag   = 1;            % [m^2] drag surface area
secondary(3).Cd       = 2.2;          % [-] shape coefficient for drag
secondary(3).A_srp    = 1;            % [m^2] SRP surface area
secondary(3).Cr       = 1.31;         % [-] shape coefficient for SRP
secondary(3).x          = [];         % [-] shape coefficient for SRP
secondary(3).covariance = [];         % [-] shape coefficient for SRP


%%
pp = struct( ...
            'orbit',     orbit, ...
            'mu',        mu, ...
            'primary',   primary, ...
            'secondary', secondary, ...
            'T',         T ...
            );
end