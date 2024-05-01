function J = pseudoHess3(dv0, G)
% DvOrder2 computes the Dv given a second order expansion of the collision
% metric.
% INPUT:  dv0 = [-] Control variables linearization point
%         G   = [-] Hihg-orer tensor of DA-coefficients of order 3
%         
% OUTPUT: J   = [-] Linearized pseudo-Hessian
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
n = length(dv0);
for j = 1:n
    J(:,j) = squeeze(dv0'*G(:,:,j));
end
end