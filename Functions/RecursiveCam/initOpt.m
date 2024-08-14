function pp = initOpt(multiple,cislunar,i)
% initOpt defines the conjunctions parameters for the polynomial 
% optimization.
% 
% INPUT:
%        multiple = [-] flag that activates multiple encounters
%        cislunar = [-] flag that activates cislunar test case and dynamics
%        i        = [-] index of the conjunction in the ESA challenge table
% 
% OUTPUT:
%        pp       = [struct] optimization paramters structure
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------

%(modifiable)
if multiple && ~cislunar
    pp   = generateInitMultiple(multiple);
%     pp   = generateInitMultipleScitech(multiple);
    pp.T         = pp.T/pp.Tsc;                                                     % [-]    (1,1) Scaled orbital period
elseif ~multiple && cislunar
    % pp = generateInitHalo();
    % pp = generateInitLyapunovJpl();
    pp = generateInitNrho('perp');
%     pp = generateInitLyapunov();
    % pp = generateInitDro();
    % pp = generateInitP3Dro();
    pp.T  = 1;                                                     % [-]    (1,1) Scaled orbital period

else
    pp   = generateInitLeo(i);
    pp.T         = pp.T/pp.Tsc;                                                     % [-]    (1,1) Scaled orbital period
end
pp.gravOrd   = 0;                                                               % Order of the gravitational model

%(not modifiable)
pp.Asc     = pp.Vsc/pp.Tsc;                                                     % [km/s2]   (1,1) Acceleration scaling constant
pp.scaling = [pp.Lsc*ones(3,1); pp.Vsc*ones(3,1)];                              % [km km/s] (6,1) Vector of scaling constants
pp.et      = 478548000/pp.Tsc;                                                  % [-]       (1,1) Scaled initial ephemeris time
pp.x_pTCA  = pp.primary.x0./pp.scaling;                                         % [-]       (6,1) Scaled primary cartesian state at TCA
pp.primary.n = pp.primary.n*pp.Tsc;
pp.secondary.n = pp.secondary.n*pp.Tsc;
for k = 1:pp.n_conj
    pp.x_sTCA(:,k) = pp.secondary(k).x0./pp.scaling;                            % [-] (6,1) Scaled kth secondary cartesian state at TCA
    pp.HBR(k)      = pp.secondary(k).HBR/pp.Lsc;                                % [-] (6,1) Scaled Hard Body Radius of kth secondary
end
for k = 1:pp.n_conj
    if cislunar; r2e_s = eye(3); r2e_p = eye(3); else                           % [-] (3,3) null rotation matrix for cislunar
        r2e_s       = rtn2eci(pp.x_sTCA(1:3,k),pp.x_sTCA(4:6,k));               % [-] (3,3) RTN to ECI rotation matrix for kth secondary in Earth Orbit 
        r2e_p       = rtn2eci(pp.x_pTCA(1:3),pp.x_pTCA(4:6));                   % [-] (3,3) RTN to ECI rotation matrix for primary in Earth Orbit 
    end
    pp.P(:,:,k) = (r2e_s*pp.secondary(k).C0(1:3,1:3)*r2e_s' + ... 
                    r2e_p*pp.primary.C0(1:3,1:3)*r2e_p')/pp.Lsc^2;              % [-] (3,3) Rotated position combined covariance matrix
    D = diag(1./pp.scaling);
    pp.Cp(:,:,k) = D*pp.primary.C0*D;
    pp.Cs(:,:,k) = D*pp.secondary(k).C0*D;
end
end