function pp = aidaInit(pp)

% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

fid = fopen('write_read/AIDA_init.dat', 'w');
fprintf(fid, '%2i\n',     pp.aidaFlag1); % atmosphere flag (1:non-rotating, 2:rotating)
fprintf(fid, '%2i\n',     pp.aidaFlag2); % SRP flag (1:no shadow, 2:Earth cylindrical shadow, 3:Earth biconical shadow, 4:Earth and Moon cylindrical shadow, 5:Earth biconical and Moon cylindrical shadow, 6:Earth and Moon biconical shadow)
fprintf(fid, '%2i\n',     pp.aidaFlag3); % third body flag (1:Moon, 2:Moon and Sun)
fprintf(fid, '%2i\n',     pp.gravOrd);
fprintf(fid, '%40.12f\n', pp.primary.mass);
fprintf(fid, '%40.12f\n', pp.primary.A_drag);
fprintf(fid, '%40.12f\n', pp.primary.Cd);
fprintf(fid, '%40.12f\n', pp.primary.A_srp);
fprintf(fid, '%40.12f\n', pp.primary.Cr);
fclose(fid);
end