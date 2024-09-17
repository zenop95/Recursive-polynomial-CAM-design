function yf = computeCtrlNlp(coeff,u,pp)
% computeDvNlp Solves the polynomial CAM optimization problem using fmincon
%
% INPUT:  coeff = [struct] Structure with coefficients of the expansion
%         u        = [-] control of the reference trajectory
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
limLo(limLo == -inf) = -1000;
y0       = reshape(u,1,m*n_man);
lb       = -ones(m*n_man,1);
ub       =  ones(m*n_man,1);
options  = optimoptions( ...
                       'fmincon',                           ...
                       'Display',               'none',     ...
                       'ConstraintTolerance',   1e-12, ...
                       'OptimalityTolerance',   1e-12, ...
                       'Algorithm',             'interior-point',      ...
                       'StepTolerance',         1e-12,      ...
                       'MaxIterations',         2e3,        ...
                       'MaxFunctionEvaluations',2e3         ...
                       );
Yf  = fmincon(@(y)minFun(y,m,n_man,pp),y0,[],[],[],[],lb,ub,@(y)polyConstr(y,coeff,limUp,limLo,pp.isEqConstr),options);
yf  = reshape(Yf,m,n_man);

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

function [c_in,c_eq] = polyConstr(y,coeff,limsUp,limsLo,isEqConstr)
    n_constr = length(coeff);
    n_eq     = sum(isEqConstr);
    n_in     = n_constr - n_eq;
    c_in     = nan(2*n_in,1);
    c_eq     = nan(n_eq,1);
    jj       = 0;
    kk       = 0;
    for j = 1:n_constr
        C       = coeff(j).C;
        E       = coeff(j).E;
        val     = 0;
        for i   = 1:length(C)
            val = val + C(i)*prod(y.^E(i,:));
        end
        if isEqConstr(j)
            jj = jj+1;
            c_eq(jj) = val - limsUp(j);
        else
            kk = kk + 1;
            c_in(kk)        = val - limsUp(j);
            c_in(kk + n_in) = val + limsLo(j);
        end

    end
end

end