function ctrl = ctrlLagrange(Deltas,DAArrays,ctrl0,k,n_man,m,n_constr)
% ctrlLagrange computes the Dv given a 4th-order expansion of the collision
% metric.
% INPUT:  Delta    = [-] Difference in the collision metric to be covered by 
%                     the Delta V 
%         DAArrays = [-] All the high-order tensors
%         ctrl0    = [-] Linearization point of the control
%         k        = [-] Order of the current truncation of the expansion
% 
% OUTPUT: ctrl     = [-] Greedy solution for the 4th-order truncations of the
%                        polynomial PoC constraint
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------
fun_handles = {@pseudoHess3,@pseudoHess4,@pseudoHess5,@pseudoHess6, ...
               @pseudoHess7,@pseudoHess8,@pseudoHess9}';                        % [-] (cell) Collect function handles for Dv up to order 9
n         = m*n_man;
grad      = nan(n_constr,n);
pseudoH   = nan(n,n,n_constr);
deltaGrad = nan(n_constr,n); 
for c = 1:n_constr
    grad(c,:)      = DAArrays{c,1};                                             % [-] (n,1) Initialize gradient with 1st-order 
    pseudoH(:,:,c) = DAArrays{c,2}*(k>1);                                       % [-] (n,n) Initialize Hessian with 2nd-order (if first-order set to zeros) 
    for j = 3:k
        pseudoH(:,:,c) = pseudoH(:,:,c) + fun_handles{j-2}(ctrl0,DAArrays{c,j});% [-] (n,n) Augmented pseudo-Hessian with j-th order contribution
    end
    deltaGrad(c,:) = ctrl0'*pseudoH(:,:,c);                                     % [-] (1,n) Pseudo-gradient contribution of orders from 2 to k
    grad(c,:)      = grad(c,:) + deltaGrad(c,:);                                % [-] (1,n) Pseudo-gradient
end
scale  = normOfVec(abs(grad'))';
Deltas = Deltas./scale;
A      = [2*eye(n), grad'; grad./scale, zeros(n_constr)];       %A(abs(A)<1e-6) = 0;
b      = [zeros(n,1); Deltas];
sol    = linsolve(A,b);
ctrl   = sol(1:end-n_constr);
end