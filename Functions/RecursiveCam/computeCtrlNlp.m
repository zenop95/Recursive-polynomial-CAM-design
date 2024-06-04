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
limUp    = pp.limUp;
limLo    = pp.limLo;
limLo(limLo == -inf) = -10;
y0       = reshape(u,1,m*n_man);
lb       = -5e-5*ones(m*n_man,1);
ub       =  5e-5*ones(m*n_man,1);
options  = optimoptions( ...
                       'fmincon',                           ...
                       'Display',               'none',     ...
                       'ConstraintTolerance',   1e-10, ...
                       'Algorithm',             'interior-point',      ...
                       'StepTolerance',         1e-12,      ...
                       'MaxIterations',         2e3,        ...
                       'MaxFunctionEvaluations',2e3         ...
                       );
Yf  = fmincon(@(y)minFun(y,m,n_man,pp),y0,[],[],[],[],lb,ub,@(y)polyConstr(y,coeff,limUp,limLo),options);
yf  = reshape(Yf,m,n_man).*scale;

function J = minFun(y,m,n_man,pp)
switch pp.objFunction
    case 'fuel'
        J = sum(normOfVec(reshape(y,m,n_man)));
    case 'energy'
        J = sum(y.^2);
    otherwise
        error('The objective function must be either fuel- or energy-optimal');
end
end

function [c_in,c_eq] = polyConstr(y,coeff,limsUp,limsLo)
    n_constr = length(coeff);
%     c_in = nan(2*n_constr,1);
    c_eq = nan(n_constr,1);
    for c = 1:n_constr
        C       = coeff(c).C;
        E       = coeff(c).E;
        val     = 0;
        for i   = 1:length(C)
            val = val + C(i)*prod(y.^E(i,:));
        end
        c_eq(c) =  val - limsUp(c);
    end
    c_in    = [];

end

end