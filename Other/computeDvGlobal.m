function yf = computeDvGlobal(lim,coeffPoC,N,m,scale)
% computeDvGlobal Solves the polynomial CAM optimization problem using a
% global optimizer
%
% INPUT:  lim      = [-] Metric limit
%         coeffPoC = [struct] Structure with coefficients of the expansion
%         N        = [-] Number of nodes for the control
%         m        = [-] Number of control variables per node
%         scale    = [-] scaling coefficients
%
% OUTPUT: Yf = [-] Optimized control vector
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
eval(['sdpvar dvSym(',num2str(m*N),',1)']);
for i = 1:N
    dvNorm(i) = sqrt(sum(dvSym(1 + m*(i-1):m*i).^2));
end

obj        = min(sum(dvNorm));                % Fuel minimization objective function
polyConstr = SymDAPolyYalmip(coeffPoC,dvSym');  % Construct polynomial from DA
% if pp.returnFlag
%     coeffDist.C(abs(coeffDist.C)<1e-4) = 0;
%     polyConstr1 = SymDAPolyYalmip(coeffDist,dvSym');  % Construct polynomial from DA
% else
%     polyConstr1 = 0;
% end
constr     = polyConstr  == lim;
%               polyConstr1 == 0];                       % Polynomial Constraint on collision metric  
%% Global Optimization
options    = sdpsettings('verbose',0,      ...
                         'warning',0,      ...
                         'solver','bmibnb', ...
                         'savesolveroutput', 0, ...
                         'debug', 0);
% options.bmibnb.pdtol   = 1e-6;
% options.bmibnb.eqtol   = 1e-6;
% options.bmibnb.vartol  = 1e-4;
% options.bmibnb.relgaptol  = 1e-2;
% options.bmibnb.absgaptol  = 1e-2;
% options.bmibnb.maxiter = 1e3;
sol = optimize(constr,obj,options);
optimize(constr,obj,options);
yf  = reshape(double(dvSym),[],N).*scale;

%% check ill-conditions in polynomial
% optimize(dvSym >= -3,min(abs(polyConstr-lim)),options);
% dv  = double(dvSym)/scaleG;

end