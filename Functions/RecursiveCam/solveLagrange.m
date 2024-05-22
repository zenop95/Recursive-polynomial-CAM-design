function [ctrl,grad] = solveLagrange(Deltas,DAArrays,ctrl0,k,n_man,m,n_constr)
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
n         = m*n_man;
grad = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);

if n_constr == 1
    gradUnit  = normalize(grad,'norm');                                         % [-] (1,n) Pseudo-gradient direction
    gradNorm  = norm(grad);                                                     % [-] (1,1) Pseudo-gradient norm
    ctrlNorm  = Deltas/gradNorm;                                                % [-] (1,1) Recursive control norm
    ctrl      = ctrlNorm*gradUnit';                                             % [-] (1,n) Recursive control solution
else    
    scale  = normOfVec(abs(grad'))';
    Deltas = Deltas./scale;
    A      = [2*eye(n), grad'; grad./scale, zeros(n_constr)];
    b      = [zeros(n,1); Deltas];
    sol    = linsolve(A,b);
    ctrl   = sol(1:end-n_constr);
end
end