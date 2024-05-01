function x = propCR3BPOde(x0,u,t,varargin)
% propCR3BPOde Uses solves the ode to propagate the
% natural motion in the CR3BP from state x0 to the final state x(t(end)).
%
% INPUT: x0   = [m][m/s] 6x1 initial state of the satellite
%        t    = [s] 1x1 final time
% 
% OUTPUT: x   = [m][m/s] 6x1 state of the final point
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
if nargin == 3 
    mu = 0.0121537; 
else 
    mu = varargin{1};
end
options = odeset('AbsTol',1e-11,'RelTol',1e-11);
x0      = [x0; u];
str     = ode78(@(t,y) cr3bp(y,mu),[0,t],x0,options);
x       = str.y(1:6,end);
end

function dy = cr3bp(y,mu)
    R = y(1:3);
    V = y(4:6);
    u = y(7:9);
    r1  = sqrt((R(1)+mu)*(R(1)+mu) + R(2)*R(2) + R(2)*R(3));
    r2  = sqrt((R(1)-(1-mu))*(R(1)-(1-mu)) + R(1)*R(1) + R(3)*R(3));
    rDot = V;
    vDot = [2*y(5) + R(1) - (1-mu)*(R(1)+mu)/r1^3 - mu*(R(1)-(1-mu))/r2^3;
           -2*y(4) + R(2) - (1-mu)*R(2)/r1^3 - mu*R(2)/r2^3;
           -(1-mu)*y(3)/r1^3 - mu*y(3)/r2^3];
    dy = [rDot; vDot + u; zeros(3,1)];
end