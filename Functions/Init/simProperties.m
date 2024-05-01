function pp = simProperties(orbit,indCase)
% simProperties initializes the simulation parameters.
%
% INPUT: orbit = [str] Defines the orbit type as LEO or GEO
%
% OUTPUT: pp = [struct] Postprocess structure
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
%% Obit generation
nxOrb = 60;
if strcmpi(orbit,'geo')
%     pp = generateInitGeo(orbit);
else
    pp = generateInitRoberto(orbit,indCase);
%     pp = generateInit(orbit);
%     pp = generateInitMultiple(orbit);
%     pp = generateInitMultipleStarlink(orbit);
%     pp = generateInitMultipleMasson(orbit);
%     pp = generateInitMultipleGmm(orbit);
%     pp = generateInitMultipleLong(orbit);
end
% pp.x2cart = @MEE2RV;
% pp.cart2x = @RV2MEE;
pp.x2cart = @null;
pp.cart2x = @null;
pp.timeSubtr = 0;
%% Most often changed parameters
n_f       = 0;                                                             % [-] (1,1) Number of orbits forward propagation from conjunction
n_b       = 1.5;                                                             % [-] (1,1) Number of orbits backward propagation from conjunction
pp.enableSmdGradConstraint = false;                                        % [bool] (1,1) Flag to toggle SMD gradient constraint
pp.flagRtn               = true;                                                   % [bool] (1,1) Reference frame in which the thrust is expressed ( 0 = ECI, 1 = RTN)
pp.gmmOrder              = 1;
pp.singleObject          = false;
pp.adaptSmdLimit         = false;
pp.toggleRefineTca       = false; 
pp.ipcConstr             = false;                                           % [bool] (1,1) Toggle risk linear constraint vs minor iterations
pp.embedLinearConstraint = false;
pp.enableSkTarget        = false;                                           % [bool] (1,1) Flag to toggle SK targeting contraint

%% Propulsive system paramters
uMax = 1e-5;                 %hi thrust                                  % [km/s^2] (1,1) Maximum control acceleration
% uMax = 1e-6;               %lo thrust                                    % [km/s^2] (1,1) Maximum control acceleration
% uMax = 2e-8;               %lo thrust                                      % [km/s^2] (1,1) Maximum control acceleration
% uMax = 5e-7;               %lo thrust                                      % [km/s^2] (1,1) Maximum control acceleration
uMin = 0*uMax;                                                             % [km/s^2] (1,1) Minimum control acceleration

%% Scaling factors and scaled initial variables
Lsc = pp.primary.a;                                                        % [km]             (1,1) Distance scaling constant
Vsc = sqrt(pp.mu/Lsc);                                                     % [km/s]           (1,1) Velocity scaling constant
Tsc = Lsc/Vsc;                                                             % [s]              (1,1) Time scaling constant
Asc = Vsc/Tsc;                                                             % [km/s^2]         (1,1) Acceleration scaling constant
pp.scaling = [Lsc*ones(3,1); Vsc*ones(3,1); Asc*ones(3,1)];                % [km;km/s;km/s^2] (9,1) Scaling vector
pp.Tsc     = Tsc;

%% Fast encounter
pp.fastEncounter         = true;                                           % [bool] (1,1) Short-term encounter flag
pp.enableBplaneAvoidance = pp.fastEncounter*true;                          % [bool] (1,1) B-plane computation flag

%% Time parameters
pp.utc      = '02/03/2015 06:00:00'; % utc at tca
pp.order  = 2;                                                             % [-] (1,1) Order of the DA propagation
pp.N_forw = round(n_f*nxOrb);                                              % [-] (1,1) Nodes after the conjunction node
pp.N_back = round(n_b*nxOrb);                                              % [-] (1,1) Nodes before the conjuntion node
pp.N      = pp.N_back + pp.N_forw + 1;                                     % [-] (1,1) Total number of nodes (with conjunction node)
pp.dt     = pp.T/nxOrb;                                                    % [s] (1,1) Time step
pp.t0     = 0;                                                             % [s] (1,1) Conjunction time where the covariance is known
pp.tb     = -pp.N_back*pp.dt;                                              % [s] (1,1) Starting time of simulation
pp.tf     = pp.t0+pp.N_forw*pp.dt;                                         % [s] (1,1) Ending time of simulation
pp.t      = linspace(pp.tb,pp.tf,pp.N)';                                   % [s] (N,1) Time nodes vector
simWindow = [-pp.N_back pp.N_forw]*pp.dt;
if pp.fastEncounter
    pp.NCA0   = pp.N_back+1; 
    pp.NCAf   = pp.N_back+1; % nodes in which CA constraint is active
else
    pp.NCA0   = 2;                                                         % [-] (1,1) Initial node in which CA constraint is active
    pp.NCAf   = pp.N;                                                      % [-] (1,1) Final node in which CA constraint is active
end
pp.NSK0   = 2;                                                             % [-] (1,1) Initial node in which SK constraint is active
pp.NSKf   = pp.N;                                                          % [-] (1,1) Final node in which SK constraint is active
etTca     = utc2et(pp.utc);
pp.et     = etTca;% - pp.dt*pp.N_back;                                     % [-] (1,1) Ephemeris time of closest approach
etStart   = etTca - pp.dt*pp.N_back;                                       % [-] (1,1) Ephemeris time of simulation start
etEnd     = etTca + pp.dt*pp.N_forw;                                       % [-] (1,1) Ephemeris time of simulation finish
utcStart  = et2utc(etStart);
utcEnd    = et2utc(etEnd);

pp.dynamics  = 'aida';
pp.nDAVars   = 9;                                                          % [-] (1,1) number of independent DA variables
pp.nDepVars  = 12;                                                         % [-] (1,1) number of dependent DA variables
pp.nu_max    = 1e-4;                                                       % [-] (1,1) Maximum allowed value for the NLI
pp.nu_m      = 1e-10;                                                       % [-] (1,1) Maximum allowed value for the NLI
pp.nu_M      = 1e-2;                                                       % [-] (1,1) Maximum allowed value for the NLI

%% Scale stuff
D         = diag(1./pp.scaling(1:6));                                      % [-] (6,6) Scaling matrix for covariance
for j = 1:length(pp.secondary)
    pp.secondary(j).HBR      = pp.secondary(j).HBR/Lsc;                    % [-] (1,1) Scaled Hard Body Radius
    pp.secondary(j).C0       = D*pp.secondary(j).C0*D;                     % [-] (6,6) Scaled initial covariance
    if isempty(pp.secondary(j).x0) && ~isempty(pp.secondary(j).relState)
        pp.secondary(j).relState = pp.secondary(j).relState./pp.scaling(1:6);  % [-] (6,1) Scaled Initial relative position
    elseif ~isempty(pp.secondary(j).x0) && isempty(pp.secondary(j).relState)
        pp.secondary(j).x0 = pp.secondary(j).x0./pp.scaling(1:6);          % [-] (6,1) Scaled Initial relative position
    end
end
pp.primary.cart0 = pp.primary.cart0./pp.scaling(1:6);                      % [-] (6,1) Scaled initial cartesian state of the secondary
pp.primary.x0 =  pp.cart2x(pp.primary.cart0);                              % [-] (6,1) Scaled initial state of the secondary
pp.T        = pp.T/Tsc;                                                    % [-] (1,1) Scaled orbital period
r2e         = rtn2eci(pp.primary.cart0(1:3),pp.primary.cart0(4:6));
R2E         = [r2e zeros(3); zeros(3) r2e];
pp.primary.C0    = R2E*D*pp.primary.C0*D*R2E';                             % [-] (6,1) Scaled initial covariance of the primary in ECI
pp.et     = pp.et/pp.Tsc;                                                  % [-] (1,1) Scaled initial ephemeris time
pp.dt     = pp.dt/pp.Tsc;                                                  % [-] (1,1) Scaled time step
pp.t      = pp.t/pp.Tsc;                                                   % [-] (N,1) Scaled time vector 
pp.uMax   = uMax/Asc;                                                      % [-] (1,1) Scaled maximum control acceleration 
pp.uMin   = uMin/Asc;                                                      % [-] (1,1) Scaled minimum control acceleration 

%% Collision risk indicator selection
pp.obj       = 'risk';                                                     % [str] Type of collision metrics ('risk','max_risk','miss_distance','2d')
pp.lim       = 1e-6;                                                       % [-] (1,1) IPC limit
% miss_distance_lim = 2;                                                                                                              %[km] miss distance limit
% pp.obj       = 'miss_distance';                                            % [str] Type of collision metrics ('risk','max_risk','miss_distance','2d')
% pp.mdCone    = true;                                                       % [str] Type of collision metrics ('risk','max_risk','miss_distance','2d')
% pp.lim = miss_distance_lim/Lsc;
% pp.PoCType   = 'Maximum';                                                % [str] Type of PoC ('Constant','Maximum',Chan','Alfano')
pp.PoCType   = 'Constant';                                                 % [str] Type of PoC ('Constant','Maximum',Chan','Alfano')
% pp.PoCType   = 'Chan';                                                   % [str] Type of PoC ('Constant','Maximum',Chan','Alfano')
% pp.PoCType   = 'Alfano';                                                 % [str] Type of PoC ('Constant','Maximum',Chan','Alfano')
pp.ipc_type   = 'Constant';                                                % [str] Type of IPC ('Constant')

%%
pp                       = switchObj(pp);

%% Optimization parameters  
pp.iterMaxMin      = 10;                                                   % [-] (1,1) Maximum number of minor iterations
pp.iterMaxMaj      = 20;                                                   % [-] (1,1) Maximum number of major iterations
pp.tolMin          = 1e-6;                                                 % [-] (1,1) Minor iterations tolerance
pp.tolMaj          = 5e-3;                                                 % [-] (1,1) Major iterations tolerance
pp.trWeight        = 0;
pp.vcWeight        = 1e5;                                                  % [-] (1,1) Weight assigned to the virtual control cone in cost function (must be orders of magnitude higher than other weights) 
pp.ctrlWeight      = 1/pp.N;%pp.uMax;                                              % [-] (1,1) Weight assigned to the control cone in cost function (set to uMax for consistency, if the uMax value is changed the total cost should stay the same)
pp.loadPreviousSol = false;                                                % [bool] (1,1) Flag to load previous solution
pp.previousSolName = 'sol';
pp.adaptTrustRegion = false;                                                      % [-] (1,1) Scaled minimum control acceleration 
if pp.ipcConstr; pp.iterMaxMin = 1; end
if pp.embedLinearConstraint == true; pp.ipcConstr = false; end

%% Toggle statation keping
pp.stationKeeping = false;                                                 % [bool] (1,1) Flag to toggle SK constraint during the maneuver
pp.altSk          = false;                                                 % [bool] (1,1) Flag to toggle SK constraint during the maneuver
pp.loadTarget     = false;                                                 % [bool] (1,1) To avoid computing the target for the same scenario multiple times
pp.skSoft         = false;                                                 % [bool] (1,1) Flag to toggle SK as soft contraint
pp.skDev          = [0.05 0.05]';                                          % [deg]  (2,1) Station keeping ll box bundaries
pp.altLim         = ([35783 35790.7])/pp.scaling(1);                       % [-]    (2,1) Station keeping altitude box bundaries
pp.skTarget       = nan(6,1);           
pp.targWeight     = 100;                                                   % [-]    (1,1) Weight assigned to the the target SK slack variables in cost function
pp.skWeight       = 100;                                                   % [-]    (1,1) Weight assigned to the slack geodetic variables in cost function for soft constraint
if strcmpi(orbit,'leo')
    pp.stationKeeping = false; 
    pp.skSoft         = false; 
end
pp.sl = pp.enableSkTarget*12;                                              % [-] (1,1) Number of slack variables to add at the tail of the optimization vector

%% Toggle smd maximum gradient constraint
pp.smdSoft                 = false;                                        % [bool] (1,1) Flag to toggle SMD gradient as a soft constraint
if ~pp.enableSmdGradConstraint; pp.smdSoft = false; end
pp.smdGradSlackWeight      = 1e-1;                                         % [-] (1,1) Weight assigned to the the SMD gradient slack variables for soft constraint 
pp.maxIpcDeviation         = 0.1;                                          % [-] (1,1) Maximum allowed IPC deviation at a distance equal to HBR from the nominal point
pp.sdmGradTrig             = 0.1;                                          % [-] (1,1) This parameter is fundamental because if chosen poorly it can lead to weird results with high DeltaV
pp.gradLimContinuation     = false;

%% Operational constraints
pp.justInTime    = false;                                                  % [bool] (1,1) Toggle just in time maneuver (no attitude control)
pp.uFix          = normalize([0;-1;0],'norm');                             % [-] (3,1) Fixed thrust direction in RTN for just-in-time CAM
pp.canAlwaysFire = true;
pp.fireWindow    = [et2utc((pp.et)*pp.Tsc);
                    et2utc((pp.et +  pp.N_back/2*pp.dt)*pp.Tsc)];

%% Toggle minimum acceleration constraint
pp.enableHomotopy     = false;                                              
pp.homotopy = struct(...
                     'eps',      1e-2, ...
                     'betaTrig', 0.1, ...
                     'u_db',     1e-5 ...
                     );
pp.improveConvergence = 0;

end
