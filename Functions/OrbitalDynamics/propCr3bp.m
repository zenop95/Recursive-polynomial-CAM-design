function x = propCr3bp(x0,u,t,varargin)
% PROPKEPODE Uses solves the ode to propagate the
% keplerian orbit from state x0 to the final state x(t(end)).
%
% INPUT: x0   = [m][m/s] 6x1 initial state of the satellite
%        t    = [s] Nx1 vector of the time nodes
% 
% OUTPUT: x   = [m][m/s] 6xN matrix of the state of the orbit at each node
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
if nargin == 4 
    mu = varargin{1}; 
else 
    mu = 0.012150668; 
end
options = odeset('AbsTol',1e-11);
x0      = [x0; u];
str     = ode45(@(t,y) grav(y,mu),[0,t],x0,options);
x       = str.y(1:6,end);

end

function dS = grav(state,mu)
    x         = state(1);
    y         = state(2);
    z         = state(3);
    V         = state(4:6);
    u         = state(7:9);
    dS(1:3,1) = V;
    dS(7:9)   = 0;

    r1  = sqrt((x+mu)^2 + y^2 + z^2);
    r2  = sqrt((x+mu-1)^2 + y^2 + z^2);
    
    dS(4) = (2*V(2) + x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3) + u(1);
    dS(5) = (-2*V(1) + y - (1-mu)*y/r1^3 - mu*y/r2^3) + u(2);
    dS(6) = (-(1-mu)*z/r1^3 - mu*z/r2^3) + u(3);
    dS(3) = 0; dS(6) = 0;
    
end