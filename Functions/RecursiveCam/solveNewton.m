function ctrl = solveNewton(Deltas,DAArrays,ctrl0,k,pp)
% solveLagrange computes the Dv at the iteration of the recursive method. It
% uses a formulation with lagrange multipliers when more than one
% constraints are used, and the greedy approach in the single constraint
% case
%
% INPUT:  Deltas    = [-] Difference in the collision metric to be covered by 
%                     the Delta V 
%         DAArrays = [-] All the high-order tensors
%         ctrl0    = [-] Linearization point of the control
%         k        = [-] Order of the current truncation of the expansion
%         n_man    = [-] Number of maneuvering nodes
%         m        = [-] Number of variables per maneuverable node
%         n_constr = [-] Number of constraints
% 
% OUTPUT: ctrl     = [-] Greedy solution for the jth-order truncations of the
%                        polynomial PoC constraint
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------
n_man    = pp.n_man;
m        = pp.m;
n        = m*n_man;

grad = DAArrays{1};
H    = DAArrays{2};
for j = 3:k
    H = H + pseudoHessian(ctrl0,DAArrays{j},j);
end

% ctrl = ctrl0 - normalize(grad,'norm')'*f/norm(grad);
err = 1;
while err>1e-8
    f  = (grad+ctrl0'*H)*ctrl0-Deltas;
    fP = grad+2*ctrl0'*H;
    fS = 2*H;
    ctrl  = ctrl0 - (f./fP)'./(1 - (f./fP.^2*fS)');
    err   = norm(ctrl-ctrl0);
    ctrl0 = ctrl;
end
end