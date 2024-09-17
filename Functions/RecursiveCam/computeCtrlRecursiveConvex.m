function [yf,iters] = computeCtrlRecursiveConvex(coeff,u,pp)
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
DAArrays = cell(n_constr,DAorder);                                              % [-] (cell) Initialize cell arrays for the DA expansion high-order tensors
DeltasUp = nan(n_constr,1);
DeltasLo = nan(n_constr,1);
for c = 1:n_constr
    constPart = coeff(c).C(all(coeff(c).E==0,2));
    if isempty(constPart); constPart = 0; end
    DeltasUp(c) = limUp(c) - constPart;                                         % [-] (1,1) Relative ditance for return residual
    DeltasLo(c) = limLo(c) - constPart;                                         % [-] (1,1) Relative ditance for return residual
    DAArrays{c,2}  = zeros(n);                                                  % [-] (cell) Initialize second-order in case only first-order is used
    for k = 1:DAorder     
        DAArrays{c,k} = buildDAArray(coeff(c).C,coeff(c).E,k);                  % [-] (cell) Build cell arrays for the DA expansion high-order tensors
    end
end

%% Solve problem with convex solver
iter = 0;                                                                       % [-] (1,1) Initialize iteration counter
iters(1) = 1;
Y0   = solveConvex(DeltasUp,DeltasLo,DAArrays,zeros(pp.m*n_man,1),1,pp,0) + u;    % [-] (n,1) 1st-order greedy solution of the polynomial constraint
for k = 2:DAorder
    err  = 1;                                                                   % [-] (1,1) Initialize convergence variable
    while err > pp.tol && iter < pp.maxIter
        iter = iter + 1;                                                        % [-] (1,1) Update iteration number
        Yp   = solveConvex(DeltasUp,DeltasLo,DAArrays,Y0,k,pp,k>2) + u;             % [-] (n,1) 1st-order greedy solution of the polynomial constraint
        err  = norm(Yp-Y0);                                                     % [-] (n,1) Compute convergence variable at iteration iter
        er(iter) = err;
        Ys(:,iter) = Yp;
        Y0 = Yp;                                                                % [-] (n,1) Update linearization point for kth-order solution
    end
    % Yord(:,k) = Y0;
    iters(k) = iter-sum(iters(1:k-1));
end

yf = reshape(Yp,[],n_man);                                               % [-] (m,N) Reshape final solution to epress it node-wise
end