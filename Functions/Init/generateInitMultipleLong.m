function pp = generateInitMultipleLong(orbit)

mu  = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = .1*diag(([4.5 10 5 0.15 0.75 0.25]/1e3).^2);  % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = 0;              % [-] eccentricity 
primary.theta0 = 0;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = 0;              % [rad] inclination
primary.n      = 1.125915e-3;    %[rad/s] mean motion
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
secondary(1).tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary(1).x0       = [];         % [-] 
secondary(1).relState = [1 0.1 0 0 -2e-3 -1e-5]';                        % [km] [km/s] Relative cartesian state at TCA
secondary(1).C0       = .9*diag(([1 13 0.5 0.15 1 0.25]/1e3).^2);          % [km^2] [km^2/s^2] Covariance at TCA
secondary(1).HBR      = primary.HBR + 0.01;         % [km]
secondary(1).mass     = 100;          % [kg] mass
secondary(1).A_drag   = 1;            % [m^2] drag surface area
secondary(1).Cd       = 2.2;          % [-] shape coefficient for drag
secondary(1).A_srp    = 1;            % [m^2] SRP surface area
secondary(1).Cr       = 1.31;         % [-] shape coefficient for SRP
secondary(1).x          = [];         % [-] shape coefficient for SRP
secondary(1).covariance = [];         % [-] shape coefficient for SRP
secondary(1).w = 1;         % [-] shape coefficient for SRP
secondary(1).cdm = true;         % [-] shape coefficient for SRP
secondary(1).ang = false;

secondary(2).tca      = 3;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary(2).x0          = [];         % [-] 
secondary(2).relState = [1 0.1 0 0 -2e-3 -1e-5]';                        % [km] [km/s] Relative cartesian state at TCA
secondary(2).C0       = 9/10*diag(([4.5 10 5 0.15 0.75 0.25]/1e3).^2);      % [km^2] [km^2/s^2] Covariance at TCA in RTN
secondary(2).HBR      = primary.HBR +  0.015;         % [km]
secondary(2).mass     = 100;          % [kg] mass
secondary(2).A_drag   = 1;            % [m^2] drag surface area
secondary(2).Cd       = 2.2;          % [-] shape coefficient for drag
secondary(2).A_srp    = 1;            % [m^2] SRP surface area
secondary(2).Cr       = 1.31;         % [-] shape coefficient for SRP
secondary(2).x          = [];         % [-] shape coefficient for SRP
secondary(2).covariance = [];         % [-] shape coefficient for SRP
secondary(2).w = 1;         
secondary(2).cdm = true;         % [-] shape coefficient for SRP
secondary(2).ang = false;

%%
pp = struct( ...
            'orbit',     orbit, ...
            'mu',        mu, ...
            'primary',   primary, ...
            'secondary', secondary, ...
            'T',         T ...
            );
end