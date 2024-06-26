function ctrl = solveConvex(DeltasUp,DeltasLo,DAArrays,ctrl0,k,pp)
% linConvexProblem Solves the linear convex optimization problem of the
% multinode trajectory optimization using the conic optimization tool of 
% MOSEK.
%
% INPUT:  Deltas   = [-] Difference in the collision metric to be covered by 
%                     the Delta V 
%         DAArrays = [-] All the high-order tensors
%         pp       = [-] Parameters structure
%
% OUTPUT: ctrl      = [-] Convex solution for the quadratic optimization
%                           recursive iteration
  
% Documentation: https://docs.mosek.com/9.2/toolbox/tutorial-cqo-shared.html 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz

% Ref: https://docs.mosek.com/latest/toolbox/tutorial-qo-shared.html
%--------------------------------------------------------------------------
%% Build pseudoH and grad
n_constr = pp.n_constr;
n_man    = pp.n_man;
m        = pp.m;
n        = m*n_man;
n_opt   = (m + 1)*n_man;

grad = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);

%% Limits
prob.bux = 1e-4*ones(n_opt,1);
prob.blx = [-1e-4*ones(n_man*m,1); zeros(n_man,1)];

%% Cost function
switch pp.objFunction
    case 'fuel'
        prob.c   = [zeros(n_man*m,1); ones(n_man,1)];                                   % [struct] Define the coefficients for the objective function
    case 'energy'
        prob.qosubi = 1:n;
        prob.qosubj = 1:n;
        prob.qoval  = ones(n,1);
    otherwise
        error('The objective function must be either fuel- or energy-optimal')
end

%% Specify linear part of constraint
prob.a   = [grad, zeros(n_constr,n_man)];

%% Specify constraint upper and lower bounds
prob.blc = DeltasLo;
prob.buc = DeltasUp;

% %% Ctrl cones
if strcmpi(pp.objFunction,'fuel')
    prob.cones.type   = zeros(1,n_man);
    prob.cones.sub    = reshape([n+1:n+n_man;(1:m)'+(0:m:n-m)],n_opt,1);
    prob.cones.subptr = 1:4:4*n_man;
end
%% Run optimization
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-11;
param.MSK_DPAR_INTPNT_TOL_PFEAS    = 1e-11;
[~,res]         = mosekopt('minimize echo(0)',prob,param);
% param.MSK_IPAR_INFEAS_REPORT_AUTO = 1;
% [~,res]         = mosekopt('minimize anapro',prob,param);
try
    minorIter.feas  = res.sol.itr.prosta;                                  % [str] Extract state of the optimization solution
    if strcmpi(minorIter.feas,'unknown')
        warning('The optimizer could not find an optimal solution')
    end
    ctrl = res.sol.itr.xx(1:n);                                      % [-] (mN,1) Extract optimized vector
%     mosekResponse(res);
catch; error(['error in MOSEK: ', res.rcodestr, '. ', res.rmsg]); 
end
end