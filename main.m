%% User-defined inputs (modifiable)
initializePath();
multiple    = 2;                                                                % [-]     (1,1) flag to activate multiple encounters test case
cislunar    = 0;                                                                % [-]     (1,1) flag to activate cislunar test case
pp          = initOpt(multiple,cislunar,1);                                     % [struc] (1,1) Initialize paramters structure with conjunction data
returnTime  = -3;                                                               % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
fireTimes   = [0.5 -0.5];                               
pp.cislunar = cislunar;
pp          = defineParams(pp,fireTimes,returnTime);                            % [-] (1,1) Include optimization paramters to parameters structure

%% Non-user defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N  = pp.N;                                                                      % [-] (1,1) Number of nodes in the propagation
n_man = pp.n_man;                                                               % [-] (1,1) Number of nodes where the maneuver can be performed
if N == 1 && pp.lowThrust; error(['The algorithm needs ' ...
        'at least two nodes to define the low-thrust window']); end

m     = 3  - 2*pp.fixedDir;  pp.m  = m;                                         % [-] (1,1)  Number of optimization variables per node
u     = zeros(m,n_man);                                                         % [-] (m,N)  Ctrl of the unperturbed trajectory

%% Propagation
tic;
[lim,coeff,timeSubtr,xBall,xRetBall] = propDA(pp.DAorder,u,0,pp);

%% Optimization
if strcmpi(pp.solvingMethod,'lagrange')
        [yF,iters] = computeCtrlActiveSet(coeff,u,pp);
        % [yF,iters,er,Ys] = computeCtrlRecursive(coeff,u,pp);
elseif strcmpi(pp.solvingMethod,'fmincon')
        yF = computeCtrlNlp(coeff,u,pp);
else
    error('The solving method should be either lagrange, or fmincon')
end

if pp.fixedDir                                                                  % [-] (3,n_man) Build control matrix node-wise in the case of fixed direction
    ctrl = yF.*pp.thrustDirections(:,1:n_man);
else
    ctrl = yF;                                                                  % [-] (3,n_man) Build control matrix node-wise in the general case
end
simTime = toc - timeSubtr;     
ctrl = pp.ctrlMax*ctrl;
%% Validation
[~,~,~,x,xRetMan,x_sec,deltaTca] = propDA(1,ctrl,1,pp);                      % Validate the solution by forward propagating and computing the real PoC
if ~pp.flagMd
    metricValPoly = 10^metricValPoly;
    lim           = 10^lim;
end                                                     
%% PostProcess
postProcess(xBall,x,x_sec,xRetMan,xRetBall,lim,ctrl,deltaTca,simTime,pp)
% plotConv