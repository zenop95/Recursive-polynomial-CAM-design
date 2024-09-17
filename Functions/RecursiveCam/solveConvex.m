function ctrl = solveConvex(DeltasUp,DeltasLo,DAArrays,ctrl0,k,pp,tr)
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
n_opt    = (m + 1)*n_man;

grad = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);

%% Limits
if tr
    for i = 1:n_man*m
        limUp(i) = min(1,ctrl0(i)  + 1e-2);
        limLo(i) = max(-1,ctrl0(i) - 1e-2);
    end
    prob.bux = [limUp'; ones(n_man,1);  ones(n_constr,1)];
    prob.blx = [limLo'; zeros(n_man,1); zeros(n_constr,1)];
else
    prob.bux = [ones(n_man*(m),1);   ones(n_man,1);     ones(n_constr,1)];
    prob.blx = [-ones(n_man*m,1);    zeros(n_man,1);    zeros(n_constr,1)];
end

% prob.bux = ones(n_opt,1);
% prob.blx = [-ones(n_man*m,1); zeros(n_man,1)];

%% Cost function
switch pp.objFunction
    case 'fuel'
        prob.c   = [zeros(n_man*m,1); ones(n_man,1); ones(n_constr,1)*1e3];     % [struct] Define the coefficients for the objective function
    case 'energy'
        prob.qosubi = 1:n;
        prob.qosubj = 1:n;
        prob.qoval  = ones(n,1);
        prob.c      = [zeros(n_man*(m+1),1); ones(n_constr,1)*1e3];             % [struct] Define the coefficients for the objective function
    otherwise
        error('The objective function must be either fuel- or energy-optimal')
end

%% Specify linear part of constraint
prob.a   = [-grad, zeros(n_constr,n_man), eye(n_constr)];

%% Specify constraint upper and lower bounds
prob.buc = -DeltasLo;
prob.blc = -DeltasUp;

% %% Ctrl cones
if strcmpi(pp.objFunction,'fuel')
    prob.cones.type   = zeros(1,n_man);
    prob.cones.sub    = reshape([n+1:n+n_man;(1:m)'+(0:m:n-m)],n_opt,1);
    prob.cones.subptr = 1:4:4*n_man;
end
%% Run optimization
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-11;
param.MSK_DPAR_INTPNT_TOL_PFEAS    = 1e-11;
% [~,res]         = mosekopt('minimize echo(0)',prob,param);
param.MSK_IPAR_INFEAS_REPORT_AUTO = 1;
[~,res]         = mosekopt('minimize anapro',prob,param);
try
    feas  = res.sol.itr.prosta;                                  % [str] Extract state of the optimization solution
    if strcmpi(feas,'PRIMAL_INFEASIBLE')
        warning('The problem is primal infeasible');
    elseif strcmpi(feas,'unknown')
        warning('The optimizer could not find an optimal solution');
    end
    ctrl = res.sol.itr.xx(1:n);                                      % [-] (mN,1) Extract optimized vector
%     mosekResponse(res);
catch; error(['error in MOSEK: ', res.rcodestr, '. ', res.rmsg]); 
end
end

% function ctrl = solveConvex(DeltasUp,DeltasLo,DAArrays,ctrl0,k,pp,tr)
% % linConvexProblem Solves the linear convex optimization problem of the
% % multinode trajectory optimization using the conic optimization tool of 
% % MOSEK.
% %
% % INPUT:  Deltas   = [-] Difference in the collision metric to be covered by 
% %                     the Delta V 
% %         DAArrays = [-] All the high-order tensors
% %         pp       = [-] Parameters structure
% %
% % OUTPUT: ctrl      = [-] Convex solution for the quadratic optimization
% %                           recursive iteration
% 
% % Documentation: https://docs.mosek.com/9.2/toolbox/tutorial-cqo-shared.html 
% %
% % Author: Zeno Pavanello, 2022
% % E-mail: zpav176@aucklanduni.ac.nz
% 
% % Ref: https://docs.mosek.com/latest/toolbox/tutorial-qo-shared.html
% %--------------------------------------------------------------------------
% %% Build pseudoH and grad
% n_constr = pp.n_constr;
% n_man    = pp.n_man;
% m        = pp.m;
% n        = m*n_man;
% n_opt    = (m + 1)*n_man;
% 
% grad = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);
% ctrls = zeros(2*n,1);
% for i = 1:n
%     if ctrl0(i) >= 0
%     ctrls(i)     = ctrl0(i);
%     else
%     ctrls(i+n) = -ctrl0(i);
%     end
% end
% 
% %% Limits
% if tr
%     for i = 1:2*n
%         limUp(i) = min(1,ctrls(i)  + 1e-4);
%         limLo(i) = min(0,ctrls(i)  - 1e-4);
%     end
%     prob.bux = [limUp'; ones(n_man,1)];
%     prob.blx = [limLo'; ones(n_man,1)];
% else
%     prob.bux = ones(n_man*(2*m+1),1);
%     prob.blx = [zeros(n_man*2*m,1); zeros(n_man,1)];
% end
% 
% % prob.bux = ones(n_opt,1);
% % prob.blx = [-ones(n_man*m,1); zeros(n_man,1)];
% 
% %% Cost function
% switch pp.objFunction
%     case 'fuel'
%         prob.c   = [zeros(n_man*2*m,1); ones(n_man,1)];                                   % [struct] Define the coefficients for the objective function
%     case 'energy'
%         prob.qosubi = 1:2*n;
%         prob.qosubj = 1:2*n;
%         prob.qoval  = ones(2*n,1);
%     otherwise
%         error('The objective function must be either fuel- or energy-optimal')
% end
% 
% %% Specify linear part of constraint
% prob.a   = [grad, -grad, zeros(n_constr,n_man);];
% 
% %% Specify constraint upper and lower bounds
% prob.blc = DeltasLo;
% prob.buc = DeltasUp;
% 
% % %% Ctrl cones
% if strcmpi(pp.objFunction,'fuel')
%     prob.cones.type   = zeros(1,n_man);
%     prob.cones.sub    = reshape([n+1:n+n_man;(1:m)'+(0:m:n-m)],n_opt,1);
%     prob.cones.subptr = 1:4:4*n_man;
% end
% %% Run optimization
% param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-11;
% param.MSK_DPAR_INTPNT_TOL_PFEAS    = 1e-11;
% [~,res]         = mosekopt('minimize echo(0)',prob,param);
% % param.MSK_IPAR_INFEAS_REPORT_AUTO = 1;
% % [~,res]         = mosekopt('minimize anapro',prob,param);
% try
%     feas  = res.sol.itr.prosta;                                  % [str] Extract state of the optimization solution
%     if strcmpi(feas,'unknown')
%         warning('The optimizer could not find an optimal solution')
%     end
%     ctrl = res.sol.itr.xx(1:n)-res.sol.itr.xx(n+1:2*n);                                      % [-] (mN,1) Extract optimized vector
% %     mosekResponse(res);
% catch; error(['error in MOSEK: ', res.rcodestr, '. ', res.rmsg]); 
% end
% end