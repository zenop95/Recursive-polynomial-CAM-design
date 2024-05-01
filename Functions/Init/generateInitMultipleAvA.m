function pp = generateInitMultipleAvA(orbit)

mu  = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
p(1).C0     = 1/2*diag(([4.5 10 5 0.15 0.75 0.25]/1e2).^2);  % [km^2] [km^2/s^2] Covariance at TCA
p(1).e      = 0;              % [-] eccentricity 
p(1).theta0 = 0;              % [rad] initial true anomaly
p(1).omega  = 0;              % [rad] argument of periapsis
p(1).RAAN   = 0;              % [rad] right ascension node longitude
p(1).inc    = 0;              % [rad] inclination
p(1).n      = 1.125915e-3;    %[rad/s] mean motion
p(1).a      = (mu/p(1).n^2)^(1/3); % [km] semimajor axis
T           = 2*pi/p(1).n;         % [s] orbital period
p(1).HBR    = 0.02;           % [km]
p(1).mass   = 500;            % [kg] mass
p(1).A_drag = 1;              % [m^2] drag surface area
p(1).Cd     = 2.2;            % [-] shape coefficient for drag
p(1).A_srp  = 1;              % [m^2] SRP surface area
p(1).Cr     = 1.31;           % [-] shape coefficient for SRP
p(1).cart0 = kepler2cartesian(p(1).a,p(1).e,p(1).RAAN,p(1).inc,p(1).omega,p(1).theta0,mu);


%% Relative initial conditions in Hill reference frame
r0      = [0.03; 0.1; 0.1];      % [km] Relative position at TCA in RTN
v0      = [0; 0; 4];         % [km/s] Relative velocity at TCA in RTN
r2e     = rtn2eci(p(1).cart0(1:3),p(1).cart0(4:6));                        % [-] (3,3) DCM RTN to ECI dor the time of conjunction
R2E     = [r2e zeros(3); zeros(3) r2e];

%% Second Primary structure
secondary = struct([]);
p(2).x0       = [p(1).cart0(1:3)+[0; 0.01; 0.1]; eulang2dcm([pi/2 0  0]','xyz')*p(1).cart0(4:6)];                                 % [-] (6,1) Secondary state in ECI at the time of conjunction
p(2).cart0    = p(2).x0;
p(2).C0       = 1/2*diag(([4.5 10 5 0.15 0.75 0.25]/1e2).^2);            % [km^2] [km^2/s^2] Covariance at TCA
p(2).HBR      = 0.012;                                                     % [km]
p(2).mass     = 200;                                                       % [kg]  mass
p(2).A_drag   = 1;                                                         % [m^2] drag surface area
p(2).Cd       = 2.2;                                                       % [-]   shape coefficient for drag
p(2).A_srp    = 1;                                                         % [m^2] SRP surface area
p(2).Cr       = 1.31;                                                      % [-]   shape coefficient for SRP

coe         = cartesian2kepler(p(2).cart0);
p(2).ecc    = coe.ecc;                                                     % [-]    eccentricity 
p(2).theta0 = coe.theta;                                                   % [rad]  initial true anomaly
p(2).omega  = coe.w;                                                       % [rad]  argument of periapsis
p(2).RAAN   = coe.RAAN;                                                    % [rad]  right ascension node longitude
p(2).inc    = coe.inc;                                                     % [rad]  inclination
p(2).n      = coe.n;                                                       %[rad/s] mean motion
p(2).a      = (mu/p(2).n^2)^(1/3);                                         % [km]   semimajor axis
%%
pp = struct( ...
            'orbit',     orbit, ...
            'mu',        mu, ...
            'primary',   p, ...
            'secondary', secondary, ...
            'T',         T ...
            );
end