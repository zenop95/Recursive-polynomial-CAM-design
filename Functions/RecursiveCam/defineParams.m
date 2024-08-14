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
pp.DAorder       = 5;                                                           % [-]   (1,1) Order of the DA polynomial expansion
pp.pocType       = 1;                                                           % [-]   (1,1) PoC type (0: Constant, 1: Chan)
% pp.objFunction   = 'fuel';
pp.objFunction   = 'energy';
% pp.solvingMethod = 'lagrange';                                                  % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'newton';                                                  % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'convex';                                                  % [str] (1,1) Optimization method (recursive, fmincon)
pp.solvingMethod = 'fmincon';                                                
pp.PoCLim        = 1e-6;                                                        % [-]   (1,1) PoC limit
% pp.PoCLim        = 1e-2;                                                        % [-]   (1,1) PoC limit
% pp.PoCLim        = (1/pp.Lsc)^2;                                               % [-]   (1,1) miss distance limit
pp.nomDist       = 0.200/pp.Lsc;                                                % [-]   (1,1) Relative distance to achieve after 1 orbit

%% Operational constraints (modifiable)
pp.flagStability    = 1; % only for Cislunar
pp.lowThrust        = 0;                                                        % [bool]   (1,1) Low-thrust flag
pp.fixedDir         = 0 + pp.flagStability*pp.cislunar;                         % [bool]   (1,1) Fixed-direction flag
pp.fixedMag         = 0;                                                        % [bool]   (1,1) Fixed-magnitude flag
pp.filterMans       = 0;                                                        % [bool]   (1,1) Filtered maneuver flag
pp.maxMagConstr     = 0;                                                        % [bool]   (1,1) Filtered maneuver flag
pp.nMans            = 5;                                                        % [bool]   (1,1) Selects how many impulses to use (only valid if filerMans==1)
thrustMagnitude     = 0.1;                                                      % [mm/s^2] (1,1) Maximum acceleration if fixedMag = true
pp.thrustMagnitude  = thrustMagnitude/pp.Asc/1e6;                               % [-]      (1,1) Scaled maximum acceleration
pp.thrustDirections = repmat([0 1 0]',1,300);                                   % [-]      (3,N) Thrust directions in RTN for consecutive impulse nodes (columnwise)
pp.flagCA           = 1;
pp.flagPoCTot       = 0;
pp.flagTanSep       = 0;
pp.flagAlt          = 0;
pp.flagReturn       = 0;
pp.flagErrReturn    = 0;
pp.flagCtrlMax      = 0;
ctrlMax             = 100;                                              % [mm/s^2 or mm/s] (1,1) Maximum acceleration/deltaV if flagCtrlMax = true
pp.ctrlMax          = ctrlMax/(pp.Asc*pp.lowThrust + pp.Vsc*~pp.lowThrust)/1e6;
% pp.ctrlMax          = 1;                                                   
%% Maneuvering times (should not be modified)
if pp.cislunar; nFire = nFire/pp.Tsc*86400; nRet = nRet/pp.Tsc*86400; end       % transform days into synodic time units
nConj      = -pp.tca_sep;                                                       % [-] (1,n_conj) Conjunction times after first TCA
pp.ns      = sort(unique([nFire, nRet, nConj]),"descend")';                     % [-] (1,N) Build time discretization
canFire    = ismember(pp.ns,nFire);                                             % [-] (1,N) 1 if the node is a firing node, 0 otherwise
isRet      = ismember(pp.ns,nRet);                                              % [-] (1,N) 1 if the node is a return node, 0 otherwise
% If low-thrust model is usde, the last node of each firing window is idle
if pp.lowThrust
    for i = 2:length(pp.ns)-1
        if abs(pp.ns(i) - pp.ns(i+1)) > abs(pp.ns(i) - pp.ns(i-1)) 
            canFire(i) = 0;
        end
    end
end
pp.N        = length(pp.ns);                                                     % [-] (1,1) Total number of nodes
pp.canFire  = canFire;                                                           % [-] (1,N) 1 if the node is a firing node, 0 otherwise
pp.isRet    = isRet;                                                             % [-] (1,N) 1 if the node is a firing node, 0 otherwise
if nRet == 0; pp.isRet = zeros(pp.N,1); end 
pp.isConj   = ismember(pp.ns,nConj);                                             % [-] (1,N) 1 if the node is a conjunction, 0 otherwise
pp.t        = pp.ns*pp.T;                                                        % [-] (1,N) Time before TCA for each node (orbits for LEO, a-dimensional time units for Cislunar)
pp.n_man    = sum(pp.canFire);                                                   % [-] (1,1) Total number of firing nodes
pp.n_constr = pp.flagCA*(1 + pp.n_conj*(pp.n_conj > 1))  + pp.flagTanSep ...
              + pp.flagAlt + 6*pp.flagReturn + pp.flagErrReturn;% + pp.n_man;

%% Reference for return
r2e_p         = rtn2eci(pp.x_pTCA(1:3),pp.x_pTCA(4:6));                         % [-] (3,3) RTN to ECI rotation matrix for primary in Earth Orbit 
% pp.xReference = pp.x_pTCA + [r2e_p*[0; 20/pp.Lsc; 0]; 0; 0; 0];
pp.xReference = pp.x_pTCA;% + [r2e_p*[pp.nomDist; 0; 0]; 0; 0; 0];

%% Targets of the constraints
limUp      = [];
limLo      = [];
if pp.flagCA
    if pp.pocType == 3
        limUp = pp.PoCLim*ones(1 + pp.n_conj*(pp.n_conj > 1),1);       
        limLo = inf(1 + pp.n_conj*(pp.n_conj > 1),1);
    else
        limUp = log10(pp.PoCLim)*ones(1 + pp.n_conj*(pp.n_conj > 1),1);       
        limLo = -inf(1 + pp.n_conj*(pp.n_conj > 1),1);
    end
end
if pp.flagTanSep;    limUp   = [limUp; -.1/pp.Lsc];     limLo = [limLo; -.2/pp.Lsc]; end
if pp.flagAlt;       limUp   = [limUp; 0];              limLo = [limLo; 0]; end
if pp.flagErrReturn; limUp   = [limUp; 0];              limLo = [limLo; 0]; end
if pp.flagReturn;    limUp   = [limUp; pp.xReference];  limLo = [limLo; pp.xReference]; end
if pp.flagCtrlMax;   limUp = [limUp; ones(pp.n_man,1)]; limLo = [limLo; zeros(pp.n_man,1)]; end 
pp.limUp = limUp;
pp.limLo = limLo;

end