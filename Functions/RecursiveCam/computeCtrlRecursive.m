function [yf,iters,er,Ys] = computeCtrlRecursive(coeff,u,pp)
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
    DeltasUp(c) = limUp(c) - constPart;                                    % [-] (1,1) Relative ditance for return residual
    DeltasLo(c) = limLo(c) - constPart;                                    % [-] (1,1) Relative ditance for return residual
    DAArrays{c,2}  = zeros(n);                                             % [-] (cell) Initialize second-order in case only first-order is used
    for k = 1:DAorder
        DAArrays{c,k} = buildDAArray(coeff(c).C,coeff(c).E,k);             % [-] (cell) Build cell arrays for the DA expansion high-order tensors
    end
end

%% Define permutations
y = [];
for j = 0:n_constr
    perm = unique(perms([ones(1,n_constr-j),zeros(1,j)]),'rows');
    y = [y; perm];
end
y(end,:)           = [];                                                        % eliminate case where all constraints are inactive
y(:,pp.isEqConstr) = 1;                                                         % always activate all equality constraints
y                  = logical(unique(y,'rows'))';
comb               = size(y,2);
J                  = nan(comb,1);
Y0s                = nan(pp.m*pp.n_man,comb);

%% Solve problem with each possible set of active constraints
for j = 1:comb
    actInd     = y(:,j);
    actCoeffs  = DAArrays(actInd,:);
    actDeltaUp = DeltasUp(actInd);     
    Y0   = solveLagrange(actDeltaUp,actCoeffs,zeros(n,1),1,actInd,pp) + u;    % [-] (n,1) 1st-order greedy solution of the polynomial constraint
    iter = 0;                                                                       % [-] (1,1) Initialize iteration counter
    iters = nan(7,1);
    iters(1) = 1;
    for k = 2:DAorder
        err  = 1;                                                               % [-] (1,1) Initialize convergence variable
        while err > pp.tol && iter < pp.maxIter
            iter = iter + 1;                                                    % [-] (1,1) Update iteration number
            Yp = solveLagrange(actDeltaUp,actCoeffs,Y0,k,actInd,pp) + u;        % [-] (n,1) kth-order Lagrange solution of the polynomial constraint
            err  = norm(Yp-Y0);                                                 % [-] (n,1) Compute convergence variable at iteration iter
            er(iter,j) = err;
            Ys(:,iter) = Yp;
            Y0 = (1-pp.alpha)*Y0 + pp.alpha*Yp;                                           % [-] (n,1) Update linearization point for kth-order solution
        end
        % Yord(:,k) = Y0;
        iters(k) = iter-sum(iters(1:k-1));
    end
    Y0s(:,j) = Y0;
    J(j) = Y0'*Y0;
    grad  = psuedoGradient(DAArrays,Y0,DAorder,n_constr,n);
    g = actInd;
    for i = 1:n_constr
        if ~actInd(i)
            g(i) = grad(i,:)*Y0 - DeltasUp(i) <= 1e-11;
        end
    end
    respConstr(j) = all(g);
end
J(~respConstr) = nan; if sum(isnan(J)) == comb; error('No feasible solution'); end
[~,c]          = min(J);
ctrlcol        = Y0s(:,c);

yf = reshape(ctrlcol,[],n_man);                                               % [-] (m,N) Reshape final solution to epress it node-wise
end