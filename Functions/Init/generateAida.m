function [pp] = generateAida(pp)
% GENERATEAIDA define here the AIDA parameters for primary and secondary satellite
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

% Generic flags and gravity order
pp.aida.flag1            = 0;                                              % atmosphere flag (1:non-rotating, 2:rotating)
pp.aida.flag2            = 0;                                              % SRP flag (1:no shadow, 2:Earth cylindrical shadow, 3:Earth biconical shadow, 4:Earth and Moon cylindrical shadow, 5:Earth biconical and Moon cylindrical shadow, 6:Earth and Moon biconical shadow)
pp.aida.flag3            = 0;                                              % third body flag (1:Moon, 2:Moon and Sun)
pp.aida.gravOrd          = 0;                                             % Order of the EGM2008 gravitational model
end

