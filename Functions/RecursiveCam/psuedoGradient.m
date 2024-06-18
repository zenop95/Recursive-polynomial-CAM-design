function [grad,pseudoH] = psuedoGradient(DAArrays,ctrl0,k,n_constr,n)
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
grad      = nan(n_constr,n);
pseudoH   = nan(n,n,n_constr);
deltaGrad = nan(n_constr,n); 
for c = 1:n_constr
    grad(c,:)      = DAArrays{c,1};                                             % [-] (n,1) Initialize gradient with 1st-order 
    pseudoH(:,:,c) = DAArrays{c,2}*(k>1);                                       % [-] (n,n) Initialize Hessian with 2nd-order (if first-order set to zeros) 
    for j = 3:k
        pseudoH(:,:,c) = pseudoH(:,:,c) + pseudoHessian(ctrl0,DAArrays{c,j},j); % [-] (n,n) Augmented pseudo-Hessian with j-th order contribution
    end
    deltaGrad(c,:) = ctrl0'*pseudoH(:,:,c);                                     % [-] (1,n) Pseudo-gradient contribution of orders from 2 to k
    grad(c,:)      = grad(c,:) + deltaGrad(c,:);                                % [-] (1,n) Pseudo-gradient
end
end