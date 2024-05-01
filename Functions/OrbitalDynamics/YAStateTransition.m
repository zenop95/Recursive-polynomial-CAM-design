function [YAMatrix] = YAStateTransition(a,e,TA0,t,t0,varargin)
% YASTATETRANSITION outputs the Yamanaka-Ankersen state transition matrix
% that can be used to compute the linearized dynamics of the relative motion
% between to spacecraft. 
% INPUT: a [m]      = semimajor axis of the orbit of the spacecraft
%        e [-]      = eccentricity of the orbit of the spacecraft
%        TA0 [rad]  = offset true anomaly
%        t [s]      = end time of the evolution of the dynamics represented by
%                     the YA state transition matrix
%        t0 [s]     = start time of the evolution of the dynamics represented by
%                     the YA state transition matrix
% OUTPUT: YKMatrix = Yamanaka-Ankersen state transition matrix
%
% Reference: Yamanaka, K., & Ankersen, F. (2002). New state transition 
%            matrix for relative motion on an arbitrary Keplerian orbit. 
%            Journal of Guidance, Control, and Dynamics, 40(11), 2917–2927.
%            https://doi.org/10.2514/1.G002723
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

if nargin > 5
    mu = varargin{1}; 
else 
    mu = 398600.4418; 
end
p       = a*(1-e^2);           %[m] semilatus rectum
h       = sqrt(mu*p);          %[m^3/s^2] specific angular momentum
n       = sqrt(mu/a^3);        %[rad/s] mean motion
tOffset = trueAnomaly2time(n,e,TA0);
TA0     = time2trueAnomaly(n,e,t0+tOffset);
TA      = time2trueAnomaly(n,e,t+tOffset);
k       = sqrt(h)/p;

% at time t0
rho0 = 1+e*cos(TA0);
s0   = rho0*sin(TA0);
c0   = rho0*cos(TA0);

% first transform is to get to the pseudo initial conditions
psTransform = [rho0*eye(3) zeros(3,3);
              -e*sin(TA0)*eye(3) 1/(k^2*rho0)*eye(3)];

% compute the pseudoinital values
psMatrix =[1-e^2 0 3*e*s0*(1/rho0+1/rho0^2) -e*s0*(1+1/rho0) 0 -e*c0+2;
           0 1-e^2 0 0 0 0;
           0 0 -3*s0*(1/rho0+(e/rho0)^2) s0*(1+1/rho0) 0 c0-2*e;
           0 0 -3*(c0/rho0+e) c0*(1+1/rho0)+e 0 -s0;
           0 0 0 0 1-e^2 0;
           0 0 3*rho0+e^2-1 -rho0^2 0 e*s0]/(1-e^2);

% compute the state at time t (passing through true anomaly change of
% variable)

% at time t
rho = (1+e*cos(TA));
s   = rho*sin(TA);
c   = rho*cos(TA);
sPr = cos(TA)+e*cos(2*TA);
cPr = -(sin(TA)+e*sin(2*TA));
J   = k^2*(t-t0);

% of difference TA-TA0
rhoD = (1+e*cos(TA-TA0));
sD   = rhoD*sin(TA-TA0);
cD   = rhoD*cos(TA-TA0);

transMatrix = [1 0 -c*(1+1/rho) s*(1+1/rho) 0 3*rho^2*J;
               0 cD/rhoD 0 0 sD/rhoD 0;
               0 0 s c 0 2-3*e*s*J;
               0 0 2*s 2*c-e 0 3*(1-2*e*s*J);
               0 -sD/rhoD 0 0 cD/rhoD 0;
               0 0 sPr cPr 0 -3*e*(sPr*J+s/rho^2)];

invPsTransform = [1/rho*eye(3) zeros(3,3); k^2*e*sin(TA)*eye(3) k^2*rho*eye(3)];

YAMatrix = invPsTransform*transMatrix*psMatrix*psTransform;

if nargin > 6 
    refFrame = lower(varargin{2});
    dcm       = [0 1 0; 0 0 -1; -1 0 0];
    Hill2Lvlh = [dcm zeros(3,3); zeros(3,3) dcm];
    switch refFrame
    case 'lvlh'
        return;
    case 'hill'
        YAMatrix = Hill2Lvlh'*YAMatrix*Hill2Lvlh;
    end
end

% Eliminate values below threshold for simplicity
YAMatrix(abs(YAMatrix)<1e-10) = 0;
end

