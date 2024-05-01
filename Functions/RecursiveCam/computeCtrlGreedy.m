function yf = computeCtrlGreedy(lim,metric,coeffPoC,u,DAorder,scale,N)
% computeCtrlGreedy Solves the polynomial CAM optimization problem using a
% recursive approach
%
% INPUT:  lim      = [-] Metric limit
%         metric   = [-] metric type
%         coeffPoC = [struct] Structure with coefficients of the expansion
%         u        = [-] control of the reference trajectory
%         DAorder  = [-] DA order
%         scale    = [-] scaling coefficients
%         N        = [-] Number of nodes for the control
%
% OUTPUT: yF       = [-] Optimized control vector
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------

u        = reshape(u,[],1);                                                     % [-] (m,n_man) expansion point for the control
n        = length(u);                                                           % [-] (1,1) Number of scalar control variables
Delta    = lim - metric;                                                        % [-] (1,1) Collision metric residual
tol      = 1e-10;                                                               % [-] (1,1) Tolerance for the successive linearizations
maxIter  = 1e3;                                                                 % [-] (1,1) Maximum number of successive linearizations
DAArraysPoC  = cell(DAorder,1);                                                 % [-] (cell) Initialize cell arrays for the DA expansion high-order tensors
for k = 1:DAorder
    DAArraysPoC{k}  =  buildDAArray(coeffPoC.C,coeffPoC.E,k);                   % [-] (cell) Build cell arrays for the DA expansion high-order tensors
end
%% First-order Dv
Y0p   = findCtrl(DAArraysPoC{1},zeros(1,n),Delta) + u;                          % [-] (n,1) 1st-order greedy solution of the polynomial constraint
Yp    = Y0p;                                                                    % [-] (n,1) Initialize linearization point for 2nd-order solution

%% Higher-orders Dv
for k = 2:DAorder
    err  = 1;                                                                   % [-] (1,1) Initialize convergence variable
    iter = 0;                                                                   % [-] (1,1) Initialize iteration counter
    while err > tol && iter < maxIter
        iter = iter + 1;                                                        % [-] (1,1) Update iteration number
        Yp  = greedyCtrl(Delta,DAArraysPoC,Y0p-u,k) + u;                        % [-] (n,1) kth-order greedy solution of the polynomial constraint
        err  = norm(Yp-Y0p);                                                    % [-] (n,1) Compute convergence variable at iteration iter
        Y0p = Yp;                                                               % [-] (n,1) Update linearization point for kth-order solution
%         e(iter) = err;
    end
end
yf = reshape(Yp,[],N).*scale;                                                   % [-] (m,N) Reshape final solution to epress it node-wise
end