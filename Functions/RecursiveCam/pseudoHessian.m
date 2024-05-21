function A = pseudoHessian(dv0, A, ord)
% pseudoHessian computes the Hessian given a high-order tensor and a Delta-v.
% INPUT:  dv0 = [-] Control variables linearization point
%         G   = [-] Hihg-orer tensor of DA-coefficients of order ord
%         ord = [-] Order of the initial tensor
%         
% OUTPUT: J   = [-] Linearized pseudo-Hessian
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
for k = ord:-1:3 
    A = tensorprod(A,dv0,k,1);
end
end