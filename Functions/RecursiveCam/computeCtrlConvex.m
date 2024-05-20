function yf = computeCtrlConvex(coeff,u,scale,pp)
% computeCtrlQuadratic Solves the polynomial CAM optimization problem using a
% recursive approach with a second-order convex formulation
%
% INPUT:  coeff    = [struct] Structure with coefficients of the expansion
%                             of the constraints
%         u        = [-] control of the reference trajectory
%         scale    = [-] scaling coefficients
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
DeltasUp   = nan(n_constr,1);
DeltasLo   = nan(n_constr,1);
for c = 1:n_constr
    DeltasUp(c) = limUp(c) - coeff(c).C(1);                                         % [-] (1,1) Relative ditance for return residual
    DeltasLo(c) = limLo(c) - coeff(c).C(1);                                         % [-] (1,1) Relative ditance for return residual
    DAArrays{c,2}  = zeros(n);                                                  % [-] (cell) Initialize second-order in case only first-order is used
    for k = 1:DAorder
        DAArrays{c,k} = buildDAArray(coeff(c).C,coeff(c).E,k);                  % [-] (cell) Build cell arrays for the DA expansion high-order tensors
    end
end
Y0p    = zeros(n,1);                                                                    % [-] (n,1) Initialize linearization point for 2nd-order solution

%% Higher-orders Dv
for k = 1:DAorder
    err  = 1;                                                                   % [-] (1,1) Initialize convergence variable
    iter = 0;                                                                   % [-] (1,1) Initialize iteration counter
    while err > tol && iter < maxIter
        iter = iter + 1;                                                        % [-] (1,1) Update iteration number
        Yp  = solveConvex(DeltasUp,DeltasLo,DAArrays,Y0p,k,pp) + u;                              % [-] (n,1) kth-order greedy solution of the polynomial constraint
        err  = norm(Yp-Y0p);                                                    % [-] (n,1) Compute convergence variable at iteration iter
        Y0p = Yp;                                                               % [-] (n,1) Update linearization point for kth-order solution
        er(iter) = err;
    end
end
yf = reshape(Yp,[],n_man).*scale;                                                   % [-] (m,N) Reshape final solution to epress it node-wise
end