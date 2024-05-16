function ctrl = findCtrl(Deltas,DAArrays,ctrl0,k,n_man,m,n_constr)
% findCtrl computes the Dv at the iteration of the recursive method. It
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