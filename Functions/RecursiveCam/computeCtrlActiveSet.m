function [yf,iters] = computeCtrlActiveSet(coeff,u,pp)
% computeCtrlRecursive Solves the polynomial CAM optimization problem using a
% recursive approach with lagrange multiplier formulation
%
% INPUT: coeff = [struct] Structure with coefficients of the expansion
%                         of the constraints
%        u     = [-]     (m,1) Nominal control (usually 0)
%        pp    = [-]     (m,1) Paramters structure
%
% OUTPUT: yf    = [-] Optimized control vector
%         iters = [-] Iterations used to reach convergence at each order

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
    DAArrays{c,2}  = zeros(n);                                                  % [cell] (n_constr,2) Initialize second-order in case only first-order is used
    for k = 1:DAorder
        DAArrays{c,k} = buildDAArray(coeff(c).C,coeff(c).E,k);                  % [cell] (n_constr,DAorder) Build cell arrays for the DA expansion high-order tensors
    end
end

%% Solve problem with active set algorithm
H            = 2*eye(n); 
f            = zeros(n,1);
grad         = psuedoGradient(DAArrays,zeros(n,1),1,n_constr,n);

% If the gradient of any constraint is null, do not use it in the first order problem (it ould make it infeasible)
indEliminate = [];
n_in         = pp.n_in;
for j = 1:size(grad,1)
    if sum(abs(grad(j,:))) == 0 
       indEliminate = [indEliminate; j];
       if ~pp.isEqConstr(j)
           n_in = n_in - 1;
       end
    end
end
DAArrays1 = DAArrays; DeltasUp1 = DeltasUp; a1 = pp.isEqConstr;
DAArrays1(indEliminate,:) = []; 
DeltasUp1(indEliminate)   = [];
a1(indEliminate)          = [];
n_constr1                 = n_constr - length(indEliminate);
grad                      = psuedoGradient(DAArrays1,zeros(n,1),k,n_constr1,n);
gradEq  = grad(a1,:);
beq     = DeltasUp1(a1); if isempty(beq); beq = zeros(0,1); end 
gradIn  = grad(~a1,:);
bin     = DeltasUp1(~a1); if isempty(bin); bin = zeros(0,1); end

% Solve first order QP
options = mpcActiveSetOptions('ConstraintTolerance',1e-11);
Y0      = mpcActiveSetSolver(H,f,gradIn,bin,gradEq,beq,false(n_in,1),options); Yp = Y0;
% Restore original constraints
beq     = DeltasUp(pp.isEqConstr); if isempty(beq); beq = zeros(0,1); end
bin     = DeltasUp(~pp.isEqConstr); 
iter = 1;
Ys(:,1) = Y0;                                                                   % [-]    (m*n_man,1) Linearization point
alpha = pp.alpha;
isConstrAct = false(pp.n_in,1);                                                 % [bool] (n_constr,1) Initialize active set
iters = nan(DAorder,1);
iters(1) = 1;
for k = 2:DAorder
    err = 1;
    while err > pp.tol && iter < pp.maxIter
        iter = iter + 1;                                                        % [-] (1,1) Update iteration number
        grad = psuedoGradient(DAArrays,Y0,k,n_constr,n);                        % [-] (n_constr,m) Define pseudo-gradient
        gradEq = grad(pp.isEqConstr,:);                                         % [-] (n_eq,m) Define equality coefficient matrix
        gradIn = grad(~pp.isEqConstr,:);                                        % [-] (n_in,m) Define inequality coefficient matrix
        grads(iter,:) = gradIn(1,:);
        for j = pp.n_in
            isConstrAct(j) = logical(abs(gradIn(j,:)*Y0 - bin(j)) < 1e-5);      % [bool] (n_constr,1) Initialize active set
        end
        Yp = mpcActiveSetSolver(H,f,gradIn,bin,gradEq,beq,isConstrAct,options); % [-] (m,1) Find solution to QP    
        err  = norm(Yp-Y0);                                                     % [-] (n,1) Compute convergence variable at iteration iter
        er(iter) = err;
        Ys(:,iter) = Yp;
        Y0 = (1-alpha)*Y0 + alpha*Yp;                                           % [-] (n,1) Update linearization point for kth-order solution
    end
    iters(k) = iter-sum(iters(1:k-1));
end
yf = reshape(Yp,[],n_man);                                                      % [-] (m,N) Reshape final solution to epress it node-wise
end