function x = propCr3bp2dNoCtrl(x0,t,varargin)
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
str     = ode45(@(t,y) grav(y,mu),[0,t],x0,options);
x       = str.y(:,end);

end

function dS = grav(state,mu)
    x         = state(1);
    y         = state(2);
    V         = state(3:4);
    dS(1:2,1) = V;

    r1  = sqrt((x+mu)^2 + y^2);
    r2  = sqrt((x+mu-1)^2 + y^2);
    
    dS(3) = (2*V(2) + x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3);
    dS(4) = (-2*V(1) + y - (1-mu)*y/r1^3 - mu*y/r2^3) ;
    
end