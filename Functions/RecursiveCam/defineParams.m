function pp = defineParams(pp,fireTimes,returnTime,pclim,mdLim,mdMan)
% defineParams defines the optimization parameters for the polynomial 
% optimization.
% 
% INPUT:
%        pp         = [struct]      Optimization parameters structure
%        fireTimes  = [-] (n_man,1) Time of firing from first TCA (in orbits of the primary)
%        returnTime = [-] (1,1)     Time of return from first TCA
% 
% OUTPUT:
%        pp         = [struct] Optimization parameters structure
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
%% Optimization parameters (modifiable)
pp.DAorder       = 5;                                                           % [-] (1,1) Order of the DA polynomial expansion
pp.pocType       = 1;                                                           % [-] (1,1) PoC type (0: Constant, 1: Chan, 2: Max)
% % pp.objFunction   = 'fuel';
pp.objFunction   = 'energy';
pp.solvingMethod = 'lagrange';                                                  % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'convex';                                                  % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'fmincon';                                                
% mdLim            = 0.5;                                                         % [km]     (1,1) Miss distance limit
pp.mdLim         = (mdLim/pp.Lsc)^2;                                            % [-]      (1,1) Non-dimensionalized miss distance limit
pp.PoCLim        = pclim;                                                        % [-]      (1,1) PoC limit for each conjunction
tol               = 1e-7;                                                       % [km/s]   (1,1) Tolerance for the successive linearizations (0.1 mm/s)
% tol               = 1e-8;                                                     % [km/s^2] (1,1) Tolerance for the successive linearizations (0.01 mm/s^2)
pp.maxIter        = 5e3;                                                        % [-]      (1,1) Maximum number of successive linearizations of the high-order tensors
pp.alpha          = 0.6;                                                        % [-]      (1,1) When linearizing the PP at iteration k, the linearization point is taken as tilde{Y} = (1-alpha)*Y_{k-2} + alpha*Y_{k-1}
%% Operational constraints (modifiable)
pp.flagMd           = mdMan;                                                        % [bool]   (1,1) Miss distance instead of PoC
pp.flagStability    = 0;                                                        % [bool]   (1,1) Fix the direction of thurst to maximize stability, only for Cislunar
pp.lowThrust        = 0;                                                        % [bool]   (1,1) Toggle low-thrust dynamics
pp.fixedDir         = 0 + pp.flagStability*pp.cislunar;                         % [bool]   (1,1) Toggle fixed-direction thrust
pp.filterMans       = 1;                                                        % [bool]   (1,1) Toggle filtering the maneuver
pp.nMans            = 1;                                                        % [bool]   (1,1) Selects how many impulses to use (only valid if filerMans==1)
pp.thrustDirections = repmat([0 1 0]',1,20);                                     % [-]      (3,N) Thrust directions in RTN for consecutive impulse nodes (columnwise)
pp.flagCA           = 1;                                                        % [bool]   (1,1) Toggle Collision Avoidance constraints
pp.flagPoCTot       = 0*~pp.flagMd;                                             % [bool]   (1,1) Toggle Total Probability of Collision constraint
pp.flagMeanSma      = 0;                                                        % [bool]   (1,1) Toggle Mean Orbital Elements Station-keeping constraints
pp.flagReturn       = 0;                                                        % [bool]   (1,1) Toggle osculating Orbital Elements Station-keeping constraint (with the exception of true anomaly)
pp.flagErrReturn    = 0;                                                        % [bool]   (1,1) Toggle Cartesian error Station-keeping constraint
ctrlMax           = 1000;                                                       % [mm/s^2 or mm/s] (1,1) Maximum acceleration/deltaV (used for scaling and for setting the limit of fmicnon's variables)
pp.ctrlMax        = ctrlMax/(pp.Asc*pp.lowThrust + pp.Vsc*~pp.lowThrust)/1e6;   % [-] (1,1) Non-dimensionalized maximum acceleration/deltaV
% pp.ctrlMax          = 1;
pp.tol = tol/pp.ctrlMax;                                                        % [km/s^2 or km/s] (1,1) Scaled Tolerance
%% Maneuvering times (should not be modified)
if pp.cislunar                                                                  % transform days into synodic time units
    fireTimes = fireTimes/pp.Tsc*86400; 
    returnTime = returnTime/pp.Tsc*86400; 
end       
nConj      = -pp.tca_sep;                                                       % [-] (1,n_conj) Conjunction times after first TCA
pp.ns      = sort(unique([fireTimes, returnTime, nConj]),"descend")';           % [-] (1,N)      Build time discretization
pp.N        = length(pp.ns);                                                    % [-] (1,1)      Total number of nodes
canFire    = ismember(pp.ns,fireTimes);                                         % [-] (1,N)      1 if the node is a firing node, 0 otherwise
isRet      = ismember(pp.ns,returnTime);                                        % [-] (1,N)      1 if the node is a return node, 0 otherwise
if pp.lowThrust                                                                 % If low-thrust model is used, the last node of each firing window is idle
    for i = 2:pp.N-1
        if abs(abs(pp.ns(i) - pp.ns(i+1)) - 0.2) > 1e-15  
            canFire(i) = 0;
        end
    end
end
pp.canFire  = canFire;                                                          % [-] (1,N) 1 if the node is a firing node, 0 otherwise
pp.canFire(end)  = 0;                                                           % [-] (1,N) Cannot fire in the last node
pp.isRet    = isRet;                                                            % [-] (1,N) 1 if the node is a return node, 0 otherwise
if returnTime == 0; pp.isRet = zeros(pp.N,1); end 
pp.isConj   = ismember(pp.ns,nConj);                                            % [-] (1,N) 1 if the node is a conjunction, 0 otherwise
pp.t        = pp.ns*pp.T;                                                       % [-] (1,N) Time before TCA for each node (orbits for LEO, non-dimensional time units for Cislunar)
pp.n_man    = sum(pp.canFire);                                                  % [-] (1,1) Total number of firing nodes

%% Build limits and define constraints
limUp         = [];                                                             % [-] (-,-) Initialize upper limits of the constraints
limLo         = [];                                                             % [-] (-,-) Initialize lwoer limits of the constraints
if pp.flagCA
    if pp.flagMd
        limUp = -pp.mdLim*ones(pp.n_conj,1);                                    % [-] (n_conj,1) miss distance upper limits; the miss distance is defined as negative to be used in the same way as PoC (constraint < 0)  
    else
        limUp = log10(pp.PoCLim)*ones(pp.flagPoCTot + ...
                                    ~pp.flagPoCTot*pp.n_conj,1);                % [-] (n_conj,1) PoC upper limits
    end
    limLo = -inf(pp.flagPoCTot + ~pp.flagPoCTot*pp.n_conj,1);                   % [-] (n_conj,1) Lower limits of CA constraints
end   
if pp.flagErrReturn
    limUp      = [limUp; 0; 0];                                                 % [-] (2+n_conj,1) Define upper limits of SK constraint
    limLo      = [limLo; 0; 0];                                                 % [-] (2+n_conj,1) Define lower limits of SK constraint
end
if pp.flagReturn
    limUp      = [limUp; zeros(5,1)];                                           % [-] (5+n_conj,1) Define upper limits of SK constraint  
    limLo = [limLo; zeros(5,1)];                                                % [-] (5+n_conj,1) Define lower limits of SK constraint
end
if pp.flagMeanSma   
    limUp = [limUp; 0; 0];                                                      % [-] (2+n_conj,1) Define upper limits of SK constraint 
    limLo = [limLo; 0; 0];                                                      % [-] (2+n_conj,1) Define lower limits of SK constraint
end
pp.limUp      = limUp;                                                          % [-] (n_constr,1) pass upper limits vector to ouptu structure      
pp.limLo      = limLo;                                                          % [-] (n_constr,1) pass lower limits vector to ouptu structure
pp.isEqConstr = logical(limUp == limLo);                                        % [-] (n_constr,1) Define if constraint is equality or inequality
pp.n_constr   = length(limUp);                                                  % [-] (1,1)        Number of polynomial constraints
pp.n_eq       = sum(pp.isEqConstr);                                                % [-] (1,1)        Number of polynomial equality constraints
pp.n_in       = pp.n_constr-pp.n_eq;                                            % [-] (1,1)        Number of polynomial inequality constraints
end