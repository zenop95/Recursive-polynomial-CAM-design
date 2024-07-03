function pp = generateInitMultipleScitech(n_conj)

mu     = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = zeros(6);       % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = 0;              % [-] eccentricity 
primary.theta0 = 0;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = deg2rad(91.673);    % [rad] inclination
primary.a      = 6800;           % [km] semimajor axis
primary.n      = (mu/primary.a^3)^(1/2);    %[rad/s] mean motion
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = 0.01;           % [km]
primary.mass   = 260;             % [kg] mass
primary.A_drag = 1;               % [m^2] drag surface area
primary.Cd     = 2.2;             % [-] shape coefficient for drag
primary.A_srp  = 1;               % [m^2] SRP surface area
primary.Cr     = 1.31;            % [-] shape coefficient for SRP
x0     = kepler2cartesian(primary.a,0,0,primary.inc,0,0,mu);

%% Secondary structure
secondary = struct();
tca_sep  = [0 3660 13200 17100]; % time between the first conjunctions and the consecutive ones (first element must be zero)
tca_sep  = tca_sep(1:n_conj); 
dx = [0.01 0.01 0 -0.01; 
      0.008 0 0.005 0.01; 
     -0.02 0.01 -0.02 0.007;
     -0.1 -0.1 0 -0.1;
     -11  0.1 -13 -10;
     0.1 11 0.1 0.1];
for j = 1:n_conj
    xBall(:,j) = propKepOde(x0,zeros(3,1),tca_sep(j)+59*60,mu);
    [r2e(:,:,j),w] = rtn2eci(xBall(1:3,j),xBall(4:6,j));
    R2E(:,:,j) = rot6(r2e(:,:,j),zeros(3,1));
end
primary.x0     = xBall(:,1);
tca_sep  = tca_sep/T; % time between the first conjunctions and the consecutive ones (first element must be zero)
covs = [2.5, 4.5, 6, 4.5;
        12,  10,  6, 7;
        10,  5,   5, 2];
for j = 1:n_conj
    % secondary(j).x0       = absState(:,j);         
    secondary(j).x0       = xBall(:,j) + R2E(:,:,j)*dx(:,j);         
    secondary(j).mass     = 260;          % [kg] mass
    secondary(j).C0       = [diag(covs(:,j)), zeros(3); zeros(3,6)];
    secondary(j).HBR      = primary.HBR + 0.005;         % [km]
    secondary(j).A_drag   = 1;            % [m^2] drag surface area
    secondary(j).Cd       = 2.2;          % [-] shape coefficient for drag
    secondary(j).A_srp    = 1;            % [m^2] SRP surface area
    secondary(j).Cr       = 1.31;         % [-] shape coefficient for SRP
end
%%

pp = struct( ...
            'mu',        mu, ...
            'Lsc',       primary.a, ...                                         % [km]   (1,1) Distance scaling constant
            'Vsc',       sqrt(mu/primary.a), ...                                % [km/s] (1,1) Velocity scaling constant
            'Tsc',       sqrt(primary.a^3/mu), ...                              % [s]    (1,1) Time scaling constant
            'primary',   primary, ...
            'secondary', secondary, ...
            'T',         T, ...
            'tca_sep',   tca_sep, ...
            'n_conj',    n_conj ...
            );
end