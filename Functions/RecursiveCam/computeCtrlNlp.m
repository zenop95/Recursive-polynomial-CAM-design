function yf = computeCtrlNlp(coeff,u,scale,pp)
% computeDvNlp Solves the polynomial CAM optimization problem using fmincon
%
% INPUT:  lims      = [-] Metric limit
%         coeffPoC = [struct] Structure with coefficients of the expansion
%         u        = [-] control of the reference trajectory
%         n_man        = [-] Number of nodes for the control
%         m        = [-] Number of control variables per node
%         scale    = [-] scaling coefficients
%
% OUTPUT: Yf = [-] Optimized control vector
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
n_man    = pp.n_man;
m        = pp.m;
lim      = pp.lim;
y0       = reshape(u,1,m*n_man);
lb       = -5e-4*ones(m*n_man,1);
ub       =  5e-4*ones(m*n_man,1);
options  = optimoptions( ...
                       'fmincon',                           ...
                       'Display',               'none',     ...
                       'ConstraintTolerance',   1e-10, ...
                       'Algorithm',             'interior-point',      ...
                       'StepTolerance',         1e-12,      ...
                       'MaxIterations',         2e3,        ...
                       'MaxFunctionEvaluations',2e3         ...
                       );
Yf  = fmincon(@(y)minFun(y,m,n_man),y0,[],[],[],[],lb,ub,@(y)polyConstr(y,coeff,lim),options);
yf  = reshape(Yf,m,n_man).*scale;

function J = minFun(y,m,n_man)
    J = sum(normOfVec(reshape(y,m,n_man)));
%     J = sum(y.^2);
end

function [c_in,c_eq] = polyConstr(y,coeff,lims)

    c_eq = nan(length(coeff),1);
    for c = 1:length(coeff)
        C       = coeff(c).C;
        E       = coeff(c).E;
        val     = 0;
        for i   = 1:length(C)
            val = val + C(i)*prod(y.^E(i,:));
        end
        c_eq(c) = val - lims(c);
    end
    c_in    = [];

end

end