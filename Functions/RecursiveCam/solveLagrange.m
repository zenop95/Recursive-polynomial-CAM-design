function ctrl = solveLagrange(DeltasUp,DeltasLo,DAArrays,ctrl0,k,pp)
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
%         n_constr = [-] Number of constraints
% 
% OUTPUT: ctrl     = [-] Greedy solution for the jth-order truncations of the
%                        polynomial PoC constraint
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------
n_constr     = pp.n_constr;
n_man        = pp.n_man;
m            = pp.m;
n            = m*n_man;
grad         = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);

%%% Check if any gradient is null at first order, in which case do not use
%%% it for the optimization
if k ==1
    indEliminate = [];
    for j = 1:size(grad,1)
        if sum(grad(j,:)) == 0 
           indEliminate = [indEliminate; j];
        end
    end
    DAArrays(indEliminate,:) = []; 
    DeltasUp(indEliminate)   = [];
    DeltasLo(indEliminate)   = [];
    n_constr                 = n_constr - length(indEliminate);
    grad                     = psuedoGradient(DAArrays,ctrl0,k,n_constr,n);
end
% normGrad = max(abs(grad),[],2);
% gradSc   = grad./normGrad;
% DeltasUp = DeltasUp./normGrad;
% DeltasLo = DeltasLo./normGrad;
gradSc = grad;
if n_constr == 1
    gradUnit  = normalize(grad,'norm');                                         % [-] (1,n) Pseudo-gradient direction
    gradNorm  = norm(grad);                                                     % [-] (1,1) Pseudo-gradient norm
    ctrlNorm  = DeltasUp/gradNorm;                                              % [-] (1,1) Recursive control norm
    ctrl      = ctrlNorm*gradUnit';                                             % [-] (1,n) Recursive control solution
else    
    if pp.equalityConstr
        A      = [2*eye(n), grad'; gradSc, zeros(n_constr)];
        b      = [zeros(n,1); DeltasUp];
        sol    = linsolve(A,b);
        ctrl   = sol(1:end-n_constr);
    
    else
        y = [];
        for j = 0:n_constr
            perm = unique(perms([ones(1,n_constr-j),zeros(1,j)]),'rows');
            y = [y; perm];
        end
        y        = boolean(y');
        y(:,end) = [];
        comb     = size(y,2);
        for j = 1:comb
            Delta    = DeltasUp(y(:,j));
            inact    = n_constr-sum(y(:,j));
            gr       = grad(y(:,j),:);
            grSc     = gradSc(y(:,j),:);
            A        = [2*eye(n), gr'; grSc, zeros(n_constr-inact)];
            b        = [zeros(n,1); Delta];
            sol      = linsolve(A,b);
            ctr(:,j) = sol(1:end-n_constr+inact);
            for kk = 1:n_constr
                g(kk)         = gradSc(kk,:)*ctr(:,j) <= DeltasUp(kk) &&  ...
                                gradSc(kk,:)*ctr(:,j) >= DeltasLo(kk);
            end
            J(j)          = ctr(:,j)'*ctr(:,j);
            respConstr(j) = all(g);
        end
        J(~respConstr) = nan;
        [~,c]          = min(J);
        ctrl           = ctr(:,c);
    end
end
end