function [ctrl,grad] = solveLagrange(DeltasUp,DAArrays,ctrl0,k,actInd,pp)
% solveLagrange computes the Dv at the iteration of the recursive method. It
% uses a formulation with lagrange multipliers when more than one
% constraints are used, and the greedy approach in the single constraint
% case
%
% INPUT:  Deltas    = [-] Difference in the collision metric to be covered by 
%                     the Delta V 
%         DAArrays = [-] All the high-order tensors
%         ctrl0    = [-] Linearization point of the control
%         j        = [-] Order of the current truncation of the expansion
%         n_constr = [-] Number of constraints
% 
% OUTPUT: ctrl     = [-] Greedy solution for the jth-order truncations of the
%                        polynomial PoC constraint
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------
n_constr     = sum(actInd);
n_man        = pp.n_man;
m            = pp.m;
n            = m*n_man;
grad         = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);
%%% Check if any gradient is null at first order, in which case do not use
%%% it for the optimization
if k == 1
    indEliminate = [];
    for j = 1:size(grad,1)
        if sum(abs(grad(j,:))) == 0 
           indEliminate = [indEliminate; j];
        end
    end
    DAArrays(indEliminate,:) = []; 
    DeltasUp(indEliminate)   = [];
    n_constr                 = n_constr - length(indEliminate);
    grad                     = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);
end
% normGrad = normOfVec(grad')';
% gradSc   = grad./normGrad;
% DeltasUp = DeltasUp./normGrad;
gradSc = grad;
if n_constr == 1
    gradUnit   = normalize(gradSc,'norm');                                      % [-] (1,n) Pseudo-gradient direction
    gradNorm   = norm(gradSc);                                                  % [-] (1,1) Pseudo-gradient norm
    ctrlNorm   = DeltasUp/gradNorm;                                             % [-] (1,1) Recursive control norm
    ctrl       = ctrlNorm*gradUnit';                                            % [-] (1,n) Recursive control solution
else    
    A        = [2*eye(n), grad'; gradSc, zeros(n_constr)];
    b        = [zeros(n,1); DeltasUp];
    sol      = linsolve(A,b);
    ctrl     = sol(1:end-n_constr);
end
end