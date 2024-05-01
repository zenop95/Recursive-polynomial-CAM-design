function [set] = generateMonteCarloInit(orbit)
%GENERATEMONTECARLOINIT Summary of this function goes here
%   Detailed explanation goes here

% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
mu  = 3.986004418e14;    % [m^3/s^2]
%Limits of the prameters
RAAN = 0;
omega = 0;
theta0 = 0;
uMax = [1e-3 1e-2];
inc = [0 pi];
if strcmpi(orbit,'geo')
    a  = [42e6 42.5e6];
    HBR  = [40 90]; %[m]
    dt = 5000;
else
    a  = [6700000 7200000];
    HBR  = [20 50]; %[m]
    dt = 250;
end
ecc  = [0 0.01];
C1   = [1 15]; C2 = [1 8]; C3 = [0.5 8];
C4 = [0.1 1]; C5 = [0.5 1.5]; C6 = [0.1 0.8];
A0 = [100 5000]; B0 = [0 1000]; y_off = [0 1000]; c = [1 15];
% generate set
HBR   = generateUniformDist(HBR,1);
uMax  = generateUniformDist(uMax,1);
inc   = generateUniformDist(inc,1);
a     = generateUniformDist(a,1);
ecc   = generateUniformDist(ecc,1);
    if a*(1-ecc)<=6378000
        ecc = 0;  %Modify orbits for which the perigee is to low
    end
C1    = generateUniformDist(C1,1);
C2    = generateUniformDist(C2,1);
C3    = generateUniformDist(C3,1);
C4    = generateUniformDist(C4,1);
C5    = generateUniformDist(C5,1);
C6    = generateUniformDist(C6,1);
A0    = generateUniformDist(A0,1);
B0    = generateUniformDist(B0,1);
y_off = generateUniformDist(y_off,1);
c     = generateUniformDist(c,1);
% set = struct('HBR',{},'uMax',{},'inc',{},'a',{},'ecc',{},'C0',{},'A0',{},...
%              'B0',{},'y_off',{},'c',{},'x_s',{},'x_d',{},'pp',{});
set.HBR = HBR; set.uMax = uMax; set.inc = inc;
set.a = a; set.ecc = ecc; set.A0 = A0;
set.B0 = B0; set.y_off = y_off; set.c = c;
set.C0 = diag([C1 C2 C3 C4 C5 C6]);
n = sqrt(mu/set.a^3); set.n = n; set.T = 2*pi/n;
A0 = set.A0;
B0 = set.B0;
alpha0  = pi/180*15*(set.c-3);
beta0   = pi/2 + alpha0;
r0 = [A0*cos(alpha0); -2*A0*sin(alpha0)+set.y_off; B0*cos(beta0)];
v0 = [-n*A0*sin(alpha0); -2*n*A0*cos(alpha0); -n*B0*sin(beta0)];   

%% Inertial state of the two objects
%inertial state of the satellite
set.x_s = kepler2cartesian(set.a,set.ecc,RAAN,set.inc,omega,theta0,mu);
r_s = set.x_s(1:3);
v_s = set.x_s(4:6);

%dcm from Hill to eci
dcm = lvlh2eci(r_s,v_s)*lvlh2hill()';
dr  = dcm*r0;
dv  = dcm*v0;
r_d = r_s - dr;
v_d = v_s - dv;
set.x_d = [r_d; v_d];

set.pp = struct( ...
                    'orbit', orbit, ...
                    'inc',   set.inc, ...
                    'ecc',   set.ecc, ...
                    'mu',    mu, ...
                    'HBR',   set.HBR, ...
                    'x_s',   set.x_s, ...
                    'x_d',   set.x_d, ...
                    'C0',    set.C0, ...
                    'A0',    set.A0, ...
                    'B0',    set.B0, ...
                    'y_off', set.y_off, ...
                    'T',     set.T, ...
                    'dt',    dt ...
                    );
set.deltaV = nan; set.valErr = nan; set.simTime = nan;
end

