function [grad] = newtonCoeffs(DAArrays,ctrl0,k)
% pseudoGradient computes the gradient given high-order tensors and a Delta-v.
% INPUT:  DAArrays = [-] All the high-order tensors
%         ctrl0    = [-] Linearization point of the control
%         k        = [-] Order of the current truncation of the expansion
%         n_constr = [-] Number of constraints
%         n        = [-] Number of independent variables        
%         
% OUTPUT: J   = [-] Linearized pseudo-gradient
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
H    = DAArrays{2}*(k>1);                                                  % [-] (n,n) Initialize Hessian with 2nd-order (if first-order set to zeros) 
grad = DAArrays{1} + 2*ctrl0'*H;                                          % [-] (1,n) Pseudo-gradient contribution of orders from 2 to k
for j = 3:k
    H    = pseudoHessian(ctrl0,DAArrays{j},j);                             % [-] (n,n) Augmented pseudo-Hessian with j-th order contribution
    grad = grad + k*ctrl0'*H;                                                     % [-] (1,n) Pseudo-gradient contribution of orders from 2 to k
end
end