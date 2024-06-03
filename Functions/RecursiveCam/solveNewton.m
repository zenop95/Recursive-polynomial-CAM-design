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
n_constr = pp.n_constr;
n_man    = pp.n_man;
m        = pp.m;
n        = m*n_man;

grad    = newtonCoeffs(DAArrays,ctrl0,k);   
gradNok = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);   

ctrl = ctrl0 - (gradNok*ctrl0 - Deltas)./grad';

end