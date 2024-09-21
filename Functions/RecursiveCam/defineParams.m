function pp = defineParams(pp,nFire,nRet)
% defineParams defines the optimization parameters for the polynomial 
% optimization.
% 
% INPUT:
%        pp = [struct] optimization paramters structure
% 
% OUTPUT:
%        pp    = [struct] optimization paramters structure
%        t_man = [-] Maneuvering times in orbit units
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
%% Optimization parameters (modifiable)
pp.DAorder       = 4;                                                           % [-] (1,1) Order of the DA polynomial expansion
pp.pocType       = 1;                                                           % [-] (1,1) PoC type (0: Constant, 1: Chan, 2: Max)
% % pp.objFunction   = 'fuel';
pp.objFunction   = 'energy';
pp.solvingMethod = 'lagrange';                                                  % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'convex';                                                  % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'fmincon';                                                
pp.mdLim         = (0.5/pp.Lsc)^2;                                                % [-] (1,1) miss distance limit
pp.PoCLim        = 1e-6;                                                        % [-] (1,1) PoC limit
pp.equalityConstr = 0;
tol               = 1e-7;                                                          % [km/s] (1,1) Tolerance for the successive linearizations (0.1 mm/s)
% tol               = 1e-9;                                                          % [km/s] (1,1) Tolerance for the successive linearizations (0.1 mm/s)
pp.maxIter        = 5e3;                                                        % [-] (1,1) Maximum number of successive linearizations
pp.alpha          = 0.6;                                                         % parameter to use previous iteration solution (0.1 when error return)
%% Operational constraints (modifiable)
pp.flagMd           = 1; % Miss distance instead of PoC
pp.flagStability    = 0; % only for Cislunar
pp.lowThrust        = 1;                                                        % [bool]   (1,1) Low-thrust flag
pp.fixedDir         = 0 + pp.flagStability*pp.cislunar;                         % [bool]   (1,1) Fixed-direction flag
pp.fixedMag         = 0;                                                        % [bool]   (1,1) Fixed-magnitude flag
pp.filterMans       = 0;                                                        % [bool]   (1,1) Filtered maneuver flag
pp.maxMagConstr     = 0;                                                        % [bool]   (1,1) Filtered maneuver flag
pp.nMans            = 5;                                                        % [bool]   (1,1) Selects how many impulses to use (only valid if filerMans==1)
thrustMagnitude     = 0.1;                                                      % [mm/s^2] (1,1) Maximum acceleration if fixedMag = true
pp.thrustMagnitude  = thrustMagnitude/pp.Asc/1e6;                               % [-]      (1,1) Scaled maximum acceleration
pp.thrustDirections = repmat([0 1 0]',1,5);                                   % [-]      (3,N) Thrust directions in RTN for consecutive impulse nodes (columnwise)
pp.flagCA           = 1;
pp.flagPoCTot       = 1*~pp.flagMd;
pp.flagMeanSma      = 1;
pp.flagReturn       = 0;
pp.flagErrReturn    = 0;
ctrlMax           = 1000;                                              % [mm/s^2 or mm/s] (1,1) Maximum acceleration/deltaV if flagCtrlMax = true
% ctrlMax           = 10;                                              % [mm/s^2 or mm/s] (1,1) Maximum acceleration/deltaV if flagCtrlMax = true
pp.ctrlMax        = ctrlMax/(pp.Asc*pp.lowThrust + pp.Vsc*~pp.lowThrust)/1e6;
% pp.ctrlMax          = 1;
pp.tol = tol/pp.ctrlMax;
%% Maneuvering times (should not be modified)
if pp.cislunar; nFire = nFire/pp.Tsc*86400; nRet = nRet/pp.Tsc*86400; end       % transform days into synodic time units
nConj      = -pp.tca_sep;                                                       % [-] (1,n_conj) Conjunction times after first TCA
pp.ns      = sort(unique([nFire, nRet, nConj]),"descend")';                     % [-] (1,N) Build time discretization
canFire    = ismember(pp.ns,nFire);                                             % [-] (1,N) 1 if the node is a firing node, 0 otherwise
isRet      = ismember(pp.ns,nRet);                                              % [-] (1,N) 1 if the node is a return node, 0 otherwise
% If low-thrust model is usde, the last node of each firing window is idle
if pp.lowThrust
    for i = 2:length(pp.ns)-1
        if abs(abs(pp.ns(i) - pp.ns(i+1)) - 0.2) > 1e-15  
            canFire(i) = 0;
        end
    end
end
pp.N        = length(pp.ns);                                                     % [-] (1,1) Total number of nodes
pp.canFire  = canFire;                                                           % [-] (1,N) 1 if the node is a firing node, 0 otherwise
pp.canFire(end)  = 0;                                                           % [-] (1,N) 1 if the node is a firing node, 0 otherwise
pp.isRet    = isRet;                                                             % [-] (1,N) 1 if the node is a firing node, 0 otherwise
if nRet == 0; pp.isRet = zeros(pp.N,1); end 
pp.isConj   = ismember(pp.ns,nConj);                                             % [-] (1,N) 1 if the node is a conjunction, 0 otherwise
pp.t        = pp.ns*pp.T;                                                        % [-] (1,N) Time before TCA for each node (orbits for LEO, a-dimensional time units for Cislunar)
pp.n_man    = sum(pp.canFire);                                                   % [-] (1,1) Total number of firing nodes

%% Targets of the constraints
limUp         = [];
limLo         = [];
isEqConstr    = [];
if pp.flagCA
    if pp.flagMd
        limUp = -pp.mdLim*ones(pp.n_conj,1);                                    % the miss distance is defined as negative to be used in the same way as PoC (constraint <0)  
    else
        limUp = log10(pp.PoCLim)*ones(pp.flagPoCTot + ~pp.flagPoCTot*pp.n_conj,1);       
    end
    limLo = -inf(pp.flagPoCTot + ~pp.flagPoCTot*pp.n_conj,1); 
    isEqConstr = zeros(pp.flagPoCTot + ~pp.flagPoCTot*pp.n_conj,1);
end   
if pp.flagErrReturn; limUp = [limUp; 0];   limLo = [limLo; 0]; isEqConstr = [isEqConstr; 1]; end
if pp.flagMeanSma;   limUp = [limUp; 0; 0];   limLo = [limLo; 0; 0]; isEqConstr = [isEqConstr; 1; 1]; end
pp.limUp      = limUp;
pp.limLo      = limLo;
pp.isEqConstr = boolean(isEqConstr);
pp.n_constr = length(limUp);
pp.n_eq = sum(isEqConstr);
pp.n_in = pp.n_constr-pp.n_eq;
end