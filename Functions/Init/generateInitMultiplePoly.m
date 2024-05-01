function pp = generateInitMultiplePoly(orbit)

mu     = 398600.4418;    % [km^3/s^2]
n_conj = 2;

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
primary.x0  = [6776.60657981063;
              -866.861692415755;
              -1150.36431998111;
               1.57704416507791;
               4.46511210712407;
               5.92540389971189];

%% Secondary structure
secondary = struct();
velAng    = deg2rad([40, 62, 12, 100, 18, 90, 94, 20, 143, 78]); % angle of rotation between primary and secondary velocity. The axis of rotation is the radial direction
angBool   = true;
absState = [6776.61802593467	-1440.51210611316	-2817.91718231172	6890.03879334383	-2817.86325949035;
-866.842252923163	4078.24903940634	3808.91684412513	-435.867557689529	-3808.91011641834;
-1150.37176534588	5412.03586588377	5054.57013291070	-578.299524009674	-5054.57785095697;
1.20808590955343	-3.48320177311947	-6.77796783561184	-0.137679677133498	6.59024315449470;
7.31434085309176	4.90313285377184	-0.556641257463550	5.17741147843255	0.106132559870787;
1.60488383482731	-4.62182858067473	-3.35916780375767	-5.54167488169749	-3.75394950021814];
cov      = load('covsStarlink').cov;
tca_sep  = [0 1817.29381948244 7747.41049358302 11573.2922188092 15590.4680302967]/T; % time between the first conjunctions and the consecutive ones (first element must be zero)
tca_sep  = tca_sep(1:n_conj);
for j = 1:n_conj
    secondary(j).x0       = absState(:,j);         
    secondary(j).mass     = 260;          % [kg] mass
    secondary(j).C0       = [cov(:,:,j) zeros(3,3); zeros(3,6)];            % [km^2] [km^2/s^2] Covariance at TCA
    secondary(j).HBR      = primary.HBR + 0.003;         % [km]
    secondary(j).A_drag   = 1;            % [m^2] drag surface area
    secondary(j).Cd       = 2.2;          % [-] shape coefficient for drag
    secondary(j).A_srp    = 1;            % [m^2] SRP surface area
    secondary(j).Cr       = 1.31;         % [-] shape coefficient for SRP
end
%%

pp = struct( ...
            'orbit',     orbit, ...
            'mu',        mu, ...
            'primary',   primary, ...
            'secondary', secondary, ...
            'T',         T, ...
            'tca_sep',   tca_sep, ...
            'n_conj',    n_conj ...
            );
end