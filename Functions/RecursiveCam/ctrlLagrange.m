function ctrl = ctrlLagrange(DeltaPoC,DeltaRet,DAArraysPoC,DAArraysRet,ctrl0,k,n_man,m)
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
gradPoC    = DAArraysPoC{1};                                                    % [-] (n,1) Initialize gradient with 1st-order 
gradRet    = DAArraysRet{1};                                                    % [-] (n,1) Initialize gradient with 1st-order 
pseudoHPoC = DAArraysPoC{2}*(k>1);                                              % [-] (n,n) Initialize Hessian with 2nd-order (if first-order set to zeros) 
pseudoHRet = DAArraysRet{2}*(k>1);                                              % [-] (n,n) Initialize Hessian with 2nd-order (if first-order set to zeros) 
for j = 3:k
    pseudoHPoC  = pseudoHPoC + fun_handles{j-2}(ctrl0,DAArraysPoC{j});          % [-] (n,n) Augmented pseudo-Hessian with j-th order contribution
    pseudoHRet  = pseudoHRet + fun_handles{j-2}(ctrl0,DAArraysRet{j});          % [-] (n,n) Augmented pseudo-Hessian with j-th order contribution
end
deltaGradPoC = ctrl0'*pseudoHPoC;                                               % [-] (1,n) Pseudo-gradient contribution of orders from 2 to k
deltaGradRet = ctrl0'*pseudoHRet;                                               % [-] (1,n) Pseudo-gradient contribution of orders from 2 to k
gradPoC      = gradPoC + deltaGradPoC;                                          % [-] (1,n) Pseudo-gradient
gradRet      = gradRet + deltaGradRet;                                          % [-] (1,n) Pseudo-gradient

n = m*n_man;
A = [2*eye(n), gradPoC', gradRet'; gradPoC, 0, 0; gradRet, 0, 0];
b = [zeros(n,1); DeltaPoC; DeltaRet];
sol = linsolve(A,b);
ctrl = sol(1:end-2);
end