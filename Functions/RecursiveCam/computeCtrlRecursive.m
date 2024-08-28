function [yf,iters] = computeCtrlRecursive(coeff,u,pp)
% computeCtrlRecursive Solves the polynomial CAM optimization problem using a
% recursive approach with lagrange multiplier formulation
%
% INPUT:  coeff    = [struct] Structure with coefficients of the expansion
%                             of the constraints
%         u        = [-] control of the reference trajectory
%         pp       = [-] Paramters structur
%
% OUTPUT: yF       = [-] Optimized control vector
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------
n_constr = pp.n_constr;
n_man    = pp.n_man;
limUp    = pp.limUp;
limLo    = pp.limLo;
DAorder  = pp.DAorder;
u        = reshape(u,[],1);                                                     % [-] (m,n_man) expansion point for the control
n        = length(u);                                                           % [-] (1,1) Number of scalar control variables
tol      = 1e-10;                                                               % [-] (1,1) Tolerance for the successive linearizations
maxIter  = 1e3;                                                                 % [-] (1,1) Maximum number of successive linearizations
DAArrays = cell(n_constr,DAorder);                                              % [-] (cell) Initialize cell arrays for the DA expansion high-order tensors
DeltasUp = nan(n_constr,1);
DeltasLo = nan(n_constr,1);
for c = 1:n_constr
    constPart = coeff(c).C(all(coeff(c).E==0,2));
    if isempty(constPart); constPart = 0; end
    DeltasUp(c) = limUp(c) - constPart;                                    % [-] (1,1) Relative ditance for return residual
    DeltasLo(c) = limLo(c) - constPart;                                    % [-] (1,1) Relative ditance for return residual
    DAArrays{c,2}  = zeros(n);                                             % [-] (cell) Initialize second-order in case only first-order is used
    for k = 1:DAorder
        DAArrays{c,k} = buildDAArray(coeff(c).C,coeff(c).E,k);             % [-] (cell) Build cell arrays for the DA expansion high-order tensors
    end
end
%% First-order Dv
switch lower(pp.solvingMethod)
    case {'lagrange','newton'}
        Y0   = solveLagrange(DeltasUp,DeltasLo,DAArrays,zeros(n,1),1,pp) + u;    % [-] (n,1) 1st-order greedy solution of the polynomial constraint
    case 'convex'
        Y0  = solveConvex(DeltasUp,DeltasLo,DAArrays,zeros(n,1),1,pp) + u;     % [-] (n,1) 1st-order convex solution of the polynomial constraint
end

%% Higher-orders Dv
iter = 0;                                                                       % [-] (1,1) Initialize iteration counter
Yord = Y0;
Yp = Y0;
iters = nan(7,1);
iters(1) = 1;
for k = 2:DAorder
    err  = 1;                                                                   % [-] (1,1) Initialize convergence variable
    DErr = 1;
    while err > tol && iter < maxIter
        iter = iter + 1;                                                        % [-] (1,1) Update iteration number
        if strcmpi(pp.solvingMethod,'lagrange')
            Yp = solveLagrange(DeltasUp,DeltasLo,DAArrays,Y0,k,pp) + u;                  % [-] (n,1) kth-order Lagrange solution of the polynomial constraint
        
        elseif strcmpi(pp.solvingMethod,'convex')
            Yp  = solveConvex(DeltasUp,DeltasLo,DAArrays,Y0,k,pp) + u;          % [-] (n,1) kth-order convex solution of the polynomial constraint
        
        elseif strcmpi(pp.solvingMethod,'newton')
            Yp = solveNewton(DeltasUp,DAArrays,Y0,k,pp) + u;                    % [-] (n,1) kth-order Lagrange solution of the polynomial constraint
        
        end

        err  = norm(Yp-Y0);                                                     % [-] (n,1) Compute convergence variable at iteration iter
        er(iter) = err;
        Ys(:,iter) = Yp;
        if iter > 1
            DErr(iter) = er(iter-1) - err;
        end
        Y0 = (1-pp.alpha)*Y0 + pp.alpha*Yp;                                           % [-] (n,1) Update linearization point for kth-order solution
    end
    Yord(:,k) = Y0;
    iters(k) = iter-sum(iters(1:k-1));
end
yf = reshape(Yp,[],n_man);                                               % [-] (m,N) Reshape final solution to epress it node-wise
end