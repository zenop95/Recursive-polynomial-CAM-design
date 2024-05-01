function yf = computeCtrlNlp(lim,coeffPoC,u,N,m,scale)
% computeDvNlp Solves the polynomial CAM optimization problem using fmincon
%
% INPUT:  lim      = [-] Metric limit
%         coeffPoC = [struct] Structure with coefficients of the expansion
%         u        = [-] control of the reference trajectory
%         N        = [-] Number of nodes for the control
%         m        = [-] Number of control variables per node
%         scale    = [-] scaling coefficients
%
% OUTPUT: Yf = [-] Optimized control vector
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

y0      = reshape(u,1,m*N);
lb      = -5e-4*ones(m*N,1);
ub      =  5e-4*ones(m*N,1);
options = optimoptions( ...
                       'fmincon',                           ...
                       'Display',               'none',     ...
                       'ConstraintTolerance',   1e-10, ...
                       'Algorithm',             'interior-point',      ...
                       'StepTolerance',         1e-12,      ...
                       'MaxIterations',         2e3,        ...
                       'MaxFunctionEvaluations',2e3         ...
                       );
Yf  = fmincon(@(y)minFun(y,m,N),y0,[],[],[],[],lb,ub,@(y)polyPoC(y,coeffPoC,lim),options);
yf  = reshape(Yf,m,N).*scale;

function J = minFun(y,m,N)
    J = sum(normOfVec(reshape(y,m,N)));
end

function [c_in,c_eq] = polyPoC(y,coeffPoC,lim)
    C       = coeffPoC.C;
    E       = coeffPoC.E;
    PoC     = 0;
    for i   = 1:length(C)
        PoC = PoC + C(i)*prod(y.^E(i,:));
    end
    c_eq    = PoC - lim;
    c_in    = PoC;
end

end