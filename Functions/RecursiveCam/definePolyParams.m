function pp = definePolyParams(pp,tMan)
% definePolyParams defines the optimization parameters for the polynomial 
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

pp.DAorder       = 2;                                                           % [-]   (1,1) Order of the DA polynomial expansion
pp.pocType       = 1;                                                           % [-]   (1,1) PoC type (0: Constant, 1: Chan)
pp.solvingMethod = 'greedy';                                                    % [str] (1,1) Optimization method (greedy,global,nlp,moment-relaxations)
% pp.solvingMethod = 'global';                                                    % [str] (1,1) Optimization method (greedy,global,nlp,moment-relaxations)
% pp.solvingMethod = 'nlp';                                                       % [str] (1,1) Optimization method (greedy,global,nlp,moment-relaxations)
% pp.solvingMethod = 'moment-relaxation';                                           % [str] (1,1) Optimization method (greedy,global,nlp,moment-relaxations)
pp.PoCLim        = 1e-6;                                                        % [-]   (1,1) PoC limit
pp.mdLim         = .3/pp.Lsc;                                                   % [-]   (1,1) Miss distance limit
pp.dyn = 0;
%% Operational constraints
pp.lowThrust        = 0;                                                        % [bool]   (1,1) Low-thrust flag
pp.fixedDir         = 0;                                                        % [bool]   (1,1) Fixed-direction flag
pp.fixedMag         = 0;                                                        % [bool]   (1,1) Fixed-magnitude flag
pp.filterMans       = 0;                                                        % [bool]   (1,1) Filtered maneuver flag
pp.nMans            = 1;                                                        % [bool]   (1,1) Selects how many impulses to use
thrustMagnitude     = 0.078;                                                     % [mm/s^2] (1,1) Maximum acceleration if fixedMag = true
pp.thrustMagnitude  = thrustMagnitude/pp.Asc/1e6;                    % [-]      (1,1) Scaled maximum acceleration
pp.thrustDirections = repmat([0 1 0]',1,300);                                    % [-]      (3,N) Thrust directions in RTN for consecutive impulse nodes (columnwise)

%% Maneuvering times
nFire      = tMan;
if pp.cislunar; nFire = nFire/4.34811305; end % transform days into synodic time units
nConj      = -pp.tca_sep;
pp.ns      = sort(unique([nFire, nConj]),"descend")';
canFire    = ismember(pp.ns,nFire);
if pp.lowThrust
    for i = 1:length(pp.ns)-1
        if abs(pp.ns(i) - pp.ns(i+1)) >= 0.4
            canFire(i) = 0;
        end
    end
end
pp.canFire = canFire;
pp.isConj  = ismember(pp.ns,nConj);
pp.t       = pp.ns*pp.T;                                                             % [-] (1,N) Maneuver time before TCA
pp.N       = length(pp.ns);
pp.n_man   = sum(pp.canFire);

end