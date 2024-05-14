function pp = defineParams(pp,tMan)
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
pp.solvingMethod = 'recursive';                                                 % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'recursiveLagrange';                                         % [str] (1,1) Optimization method (recursive, fmincon)
% pp.solvingMethod = 'fmincon';                                                
pp.PoCLim        = 1e-6;                                                        % [-]   (1,1) PoC limit

%% Operational constraints (modifiable)
pp.lowThrust        = 0;                                                        % [bool]   (1,1) Low-thrust flag
pp.fixedDir         = 0;                                                        % [bool]   (1,1) Fixed-direction flag
pp.fixedMag         = 0;                                                        % [bool]   (1,1) Fixed-magnitude flag
pp.filterMans       = 0;                                                        % [bool]   (1,1) Filtered maneuver flag
pp.nMans            = 1;                                                        % [bool]   (1,1) Selects how many impulses to use (only valid if filerMans==1)
thrustMagnitude     = 0.026;                                                    % [mm/s^2] (1,1) Maximum acceleration if fixedMag = true
pp.thrustMagnitude  = thrustMagnitude/pp.Asc/1e6;                               % [-]      (1,1) Scaled maximum acceleration
pp.thrustDirections = repmat([0 1 0]',1,300);                                   % [-]      (3,N) Thrust directions in RTN for consecutive impulse nodes (columnwise)

%% Maneuvering times (should not be modified)
nFire      = tMan;
if pp.cislunar; nFire = nFire/4.34811305; end                                   % transform days into synodic time units
nConj      = -pp.tca_sep;                                                       % [-] (1,n_conj) conjunction times after first TCA
pp.ns      = sort(unique([nFire, nConj]),"descend")';                           % [-] (1,N) Build time discretization
canFire    = ismember(pp.ns,nFire);                                             % [-] (1,N) 1 if the node is a firing node, 0 otherwise
% If low-thrust model is usde, the last node of each firing window is idle
if pp.lowThrust
    for i = 2:length(pp.ns)-1
        if abs(pp.ns(i) - pp.ns(i+1)) > abs(pp.ns(i) - pp.ns(i-1)) 
            canFire(i) = 0;
        end
    end
end
pp.N       = length(pp.ns);                                                     % [-] (1,1) Total number of nodes
pp.canFire = canFire;                                                           % [-] (1,N) 1 if the node is a firing node, 0 otherwise
pp.isConj  = ismember(pp.ns,nConj);                                             % [-] (1,N) 1 if the node is a conjunction, 0 otherwise
pp.t       = pp.ns*pp.T;                                                        % [-] (1,N) Time before TCA for each node (orbits for LEO, a-dimensional time units for Cislunar)
pp.n_man   = sum(pp.canFire);                                                   % [-] (1,1) Total number of firing nodes

end