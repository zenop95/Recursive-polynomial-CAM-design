function polyPoC = SymDAPolyYalmip(coeff,y)
% SymDAPolyYalmip Builds PoC symbolic polynomial for Yalmip's global optimizer
%
% INPUT:  coeff = [struct] Structure with coefficients of the expansion
%         y     = [-] Vecotr of symbolic control variables
%
% OUTPUT: polyPoC = [sym] Polynomial of the PoC
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
C = coeff.C;
E = coeff.E;
sdpvar polyPoC;
polyPoC = 0*y(1);
for i = 1:length(C)
    polyPoC = polyPoC + C(i)*prod(y.^E(i,:));
end
end