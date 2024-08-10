function ctrl = solveLagrange(Deltas,DAArrays,ctrl0,k,pp)
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
n_constr = pp.n_constr;% - (k==1)*pp.n_man;
% Deltas   = Deltas(1:n_constr);
n_man    = pp.n_man;
m        = pp.m;
n        = m*n_man;
grad     = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);
if n_constr == 1
    gradUnit  = normalize(grad,'norm');                                         % [-] (1,n) Pseudo-gradient direction
    gradNorm  = norm(grad);                                                     % [-] (1,1) Pseudo-gradient norm
    ctrlNorm  = Deltas/gradNorm;                                                % [-] (1,1) Recursive control norm
    ctrl      = ctrlNorm*gradUnit';                                             % [-] (1,n) Recursive control solution
else    
%     scale  = [1; 1e-15];%normOfVec(abs(grad'))';
    Deltas = Deltas;
    A      = [2*eye(n), grad'; grad, zeros(n_constr)];
    b      = [zeros(n,1); Deltas];
    sol    = linsolve(A,b);
    ctrl   = sol(1:end-n_constr);
%     y = [];
%     for j = 0:n_constr
%         perm = unique(perms([ones(1,n_constr-j),zeros(1,j)]),'rows');
%         y = [y; perm];
%     end
%     y        = boolean(y');
%     y(:,end) = [];
%     comb     = size(y,2);
%     for j = 1:comb
%         Delta  = Deltas(y(:,j));
%         inact  = n_constr-sum(y(:,j));
%         gr     = grad(y(:,j),:);
%         A      = [2*eye(n), gr'; gr, zeros(n_constr-inact)];
%         b      = [zeros(n,1); Delta];
%         sol    = linsolve(A,b);
%         ctr    = sol(1:end-n_constr+inact);
%         J(j)   = ctr'*ctr;
%         if j > 1
%             for kk = 1:n_constr
%                 g(kk) = grad(kk,:)*ctr <= Deltas(kk);
%             end
%             if J(j) < J(j-1) && all(g)
%                 ctrl = ctr;
%             end
%         else     
%             ctrl = ctr;
%         end
%     end
% end
end