function pp = generateInitMultipleMasson(orbit)

mu  = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = 0*diag(([0.01 0.1 0.1 0.0005 0.0075 0.001]/1e3).^2);  % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = 0.00145;              % [-] eccentricity 
primary.theta0 = pi/2;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = deg2rad(86.4);            % [rad] inclination
primary.a      = 7158;    %[rad/s] mean motion
primary.n      = (mu/primary.a^3)^(1/2); % [km] semimajor axis
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = 0.02;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0  = kepler2cartesian(primary.a,primary.e,primary.RAAN,primary.inc,primary.omega,primary.theta0,mu);

%% Secondary structure
secondary = struct();
relState = [-0.0079 0.035 0.08 0.0057 12.5526 5.489; 
            -0.0168 0.0018 0.041 -0.038 -14.9108 -0.7855; 
            0.009 -0.0352 0.003 0.006 -0.005 0.0058 
            ]';
cov      = [4.5 10 5 0.15 0.75 0.25;
            2.5 12 10 0.85 0 0.25;
            6 6 5 0.15 0.75 0.25].^2';
tca      = [27; 54; 234];
for j = 1:3
    secondary(j).tca      = tca(j);                                      % [s] TCA of conjunction w.r.t. initial time t0 = 0
    secondary(j).x0       = [];         
    secondary(j).mass     = 100;          % [kg] mass
    secondary(j).relState = relState(:,j);                                 % [km] [km/s] Relative cartesian state at TCA
    secondary(j).C0       = diag((cov(:,j)/1e3).^2);      % [km^2] [km^2/s^2] Covariance at TCA
    secondary(j).HBR      = primary.HBR + 0.01;         % [km]
    secondary(j).A_drag   = 1;            % [m^2] drag surface area
    secondary(j).Cd       = 2.2;          % [-] shape coefficient for drag
    secondary(j).A_srp    = 1;            % [m^2] SRP surface area
    secondary(j).Cr       = 1.31;         % [-] shape coefficient for SRP
    secondary(j).x          = [];         % [-] shape coefficient for SRP
    secondary(j).covariance = [];         % [-] 
    secondary(j).w = 1;        
    secondary(j).cdm = true;        
end
% save('sec','secondary')
% secondary = load('sec').secondary;  
%%

pp = struct( ...
            'orbit',     orbit, ...
            'mu',        mu, ...
            'primary',   primary, ...
            'secondary', secondary, ...
            'T',         T ...
            );
end