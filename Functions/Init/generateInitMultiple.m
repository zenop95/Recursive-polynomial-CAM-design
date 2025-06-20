function pp = generateInitMultiple(n_conj)

mu     = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = zeros(6);       % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = 0;              % [-] eccentricity 
primary.theta0 = 0;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = deg2rad(53);    % [rad] inclination
primary.a      = 6928;           % [km] semimajor axis
primary.n      = (mu/primary.a^3)^(1/2);    %[rad/s] mean motion
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = 0.003;           % [km]
primary.mass   = 260;             % [kg] mass
primary.A_drag = 1;               % [m^2] drag surface area
primary.Cd     = 2.2;             % [-] shape coefficient for drag
primary.A_srp  = 1;               % [m^2] SRP surface area
primary.Cr     = 1.31;            % [-] shape coefficient for SRP
primary.x0     = [6776.60657981063;
                 -866.861692415755;
                 -1150.36431998111;
                  1.57704416507791;
                  4.46511210712407;
                  5.92540389971189];
primary.n = 1;
%% Secondary structure
load("./data/dataConjunctionsESA.mat");
B   = table2array(data); clear data;

secondary = struct();
absState = [6776.528    -1440.312	-2817.917	 6890.038	-2817.863;
           -866.8422	 4078.39	 3808.916	-435.8675	-3808.910;
           -1150.571	 5412.035	 5054.62	-578.2995	-5054.577;
           -1.208085	  3.483201	6.777967	-0.137679	 6.590243;
           -5.314340	 -4.903132	-0.556641	 5.177411	 0.1061327;
            4.604883	  4.621828	3.359167	-5.541674   -3.753949];
xBall = [6776.60657980054	-1440.41220335898	-2817.87147274393	6890.04769213913	-2817.87145279969;
        -866.861692322414	4078.26364332834	3808.91311518611	-435.818300995367	-3808.91312053238;
        -1150.36431985724	5412.03864908012	5054.59842550835	-578.350419503152	-5054.59843260308;
        1.57704416491147	-7.41940950754972	-6.92939168524051	0.792865476269878	6.92939169495346;
        4.46511210715953	-0.949088882882253	-1.85669802153613	4.53985855706779	-1.85669800839242;
        5.92540389975895	-1.25948348728602	-2.46392149479167	6.02459578904420	-2.46392147734939];
tca_sep  = [0 1817.29381948244 7747.41049358302 11573.2922188092 15590.4680302967]/T; % time between the first conjunctions and the consecutive ones (first element must be zero)
conj_id = [1 4 5];
tca_sep = tca_sep(conj_id(1:n_conj));
Cs(:,:,1) = diag([1e-3 0.8 4e-5]);
Cs(:,:,5) = diag([0.01 0.05 5e-4]);
Cs(:,:,4) = diag([2e-3 0.6 4e-4]);
for j = 1:n_conj
    secondary(j).x0       = absState(:,conj_id(j));         
    secondary(j).mass     = 260;          % [kg] mass
    secondary(j).C0       = [Cs(:,:,conj_id(j)) zeros(3,3); zeros(3,6)];
    secondary(j).HBR      = primary.HBR + 0.003;         % [km]
    secondary(j).A_drag   = 1;            % [m^2] drag surface area
    secondary(j).Cd       = 2.2;          % [-] shape coefficient for drag
    secondary(j).A_srp    = 1;            % [m^2] SRP surface area
    secondary(j).Cr       = 1.31;         % [-] shape coefficient for SRP
    secondary(j).n = 1;
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