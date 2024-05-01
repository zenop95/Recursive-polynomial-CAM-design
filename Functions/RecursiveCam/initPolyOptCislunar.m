function pp = initPolyOptCislunar(conj)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
pp         = generateInitCislunar(conj);
pp         = generateAida(pp);                                                  % [km/s^2]         (1,1) Acceleration scaling constant
pp.mu      = 403504.497;
pp.Lsc     = 384405;                                                      % [km]             (1,1) Distance scaling constant
% pp.Tsc     = 2*pi*sqrt(pp.Lsc^3/(pp.mu)) ; % [s]              (1,1) Time scaling constant
pp.Tsc     = 375677; % [s]              (1,1) Time scaling constant
pp.Vsc     = pp.Lsc/pp.Tsc;                                                % [km/s]           (1,1) Velocity scaling constant
pp.Asc     = pp.Vsc/pp.Tsc; 
pp.scaling = [pp.Lsc*ones(3,1); pp.Vsc*ones(3,1)];
pp.et      = 478548000/pp.Tsc;                                                  % [-] (1,1) Scaled initial ephemeris time
pp.x_pTCA  = pp.primary.x0./pp.scaling;
pp.x_sTCA  = pp.secondary.x0./pp.scaling;
pp.HBR     = pp.secondary.HBR/pp.Lsc;
pp.P       = (pp.primary.C0(1:3,1:3) + pp.secondary.C0(1:3,1:3))/pp.Lsc^2;
pp.maxIter = 1e3;
end