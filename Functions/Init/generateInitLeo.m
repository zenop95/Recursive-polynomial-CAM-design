function pp = generateInitLeo(ind)

mu     = 398600.4418;    % [m^3/s^2]
load("./data/dataConjunctionsESA.mat");
B = table2array(data); clear data;

%% Primary
x0p = toColumn(B(ind,3:8));
primary = cartesian2kepler(x0p,mu);
primary.x0 = x0p;
primary.C0 = [[B(ind,9) B(ind,12)  B(ind,13);
           B(ind,12)  B(ind, 10) B(ind,14);
           B(ind,13)  B(ind, 14)  B(ind,11)] zeros(3,3); 
                                            zeros(3,6)];
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = B(ind,2)/2;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP

%% Secondary
x0s                = toColumn(B(ind,15:20)); % [km] [km/s] Secondary initial state in ECI
secondary.x0       = x0s;         
r2e = rtn2eci(x0p(1:3),x0p(4:6));
secondary.relState = [r2e' zeros(3); zeros(3) r2e']*(x0s-x0p);                          % [km] [km/s] Relative cartesian state at TCA
secondary.C0    = [[B(ind,21) B(ind,24) B(ind,25);
           B(ind,24)  B(ind,22) B(ind,26);
           B(ind,25)  B(ind,26) B(ind,23)] zeros(3,3); 
                                            zeros(3,6)];
a = cartesian2kepler(x0s);
secondary.n = a.n;
secondary.HBR        = B(ind,2)/2 + primary.HBR;         % [km]
secondary.mass       = 100;          % [kg] mass
secondary.A_drag     = 1;            % [m^2] drag surface area
secondary.Cd         = 2.2;          % [-] shape coefficient for drag
secondary.A_srp      = 1;            % [m^2] SRP surface area
secondary.Cr         = 1.31;         % [-] shape coefficient for SRP
%% Relative initial covariance in ECI reference frame

pp = struct( ...
            'mu',        mu, ...
            'Lsc',       primary.a, ...                                         % [km]   (1,1) Distance scaling constant
            'Vsc',       sqrt(mu/primary.a), ...                                % [km/s] (1,1) Velocity scaling constant
            'Tsc',       sqrt(primary.a^3/mu), ...                              % [s]    (1,1) Time scaling constant
            'primary',   primary, ...
            'secondary', secondary, ...
            'n_conj',    1, ...
            'tca_sep',   0, ...
            'T',         T ...
            );
end