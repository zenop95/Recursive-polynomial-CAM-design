beep off
format longG
close all
clear
addpath(genpath('.\data'))
addpath(genpath('.\Functions'))
addpath(genpath('.\CppExec'))
addpath(genpath('.\Path'));
addpath(genpath('C:\Program Files\Mosek\10.0\toolbox\r2017a'));
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultUicontrolFontName','Times', 'DefaultUicontrolFontSize', 14);
set(0,'DefaultUitableFontName','Times', 'DefaultUitableFontSize', 14);
set(0,'DefaultTextFontName','Times', 'DefaultTextFontSize', 14);
set(0,'DefaultUipanelFontName','Times', 'DefaultUipanelFontSize', 14);
set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultfigurecolor',[1 1 1])

%% Initialization variables
multiple = 0;                                                                   % [-]     (1,1) flag to activate multiple encounters test case
cislunar = 0;                                                                   % [-]     (1,1) flag to activate cislunar test case
pp = initPolyOpt(multiple,cislunar,1);                                          % [struc] (1,1) Initialize paramters structure with conjunction data
pp.cislunar = cislunar;
nFire = 2.5;
% if pp.cislunar
%     nFire = 3.6*86400/pp.Tsc;
% end
% nFire      = linspace(0.49,0.51,2);
% nFire      = [nFire, nFire + 1, nFire + 2];
pp = definePolyParams(pp,nFire);                                                % [-] (1,1) Include optimization paramters to parameters structure
N  = pp.N;                                                                      % [-] (1,1) Number of nodes in the propagation
n_man = pp.n_man;                                                               % [-] (1,1) Number of nodes where the maneuver can be performed
if N == 1 && pp.lowThrust; error(['The algorithm needs ' ...
        'at least two nodes to define the low-thrust window']); end
if pp.fixedMag && pp.fixedDir; error(['Both magnitude and direction ' ...
        '                  cannot be fixed at the same time']); end

m     = 3 - pp.fixedMag - 2*pp.fixedDir;                                        % [-] (1,1)  Number of optimization variables per node
pp.m  = m;
u     = zeros(m,n_man);                                                         % [-] (m,N)  Ctrl of the unperturbed trajectory
scale = ones(m,n_man);                                                          % [-] (m,N)  ~1 if polynomial scaling is used
ctrl  = nan(3,n_man);                                                           % [-] (3,N)  Initialized ctrl of the optimized trajectory

%% Propagation
% First propagation to compute convergence radius of DA variables
% if pp.DAorder > 1 && strcmpi(pp.solvingMethod,'global') % convergence radius is not defined for order 1  
%     [~,~,~,~,~,~,~,~,scale] = propDAPoly(pp.DAorder,u,scale,0,1,pp);
% end
% scale = min(min(scale))*ones(m,n_man);
%     scale = scale*1e-3;
% Second propagation to compute the DA maps for the collision metric
timeSubtr1 = 0;
tic
if pp.filterMans
    [~,~,coeffPoC,timeSubtr1] = propDAPoly(1,u,scale,0,0,pp);
    gradVec = buildDAArray(coeffPoC.C,coeffPoC.E,1);
    for j = 1:n_man
        grads(j) = norm(gradVec(1+m*(j-1):m*j));
    end
    [~,thrustNode] = sort(grads,'descend');
    thrustNode = thrustNode(1:pp.nMans);
    nFire      = pp.ns(thrustNode)';
    nConj     = -pp.tca_sep;
    pp.ns      = sort(unique([nFire, nConj]),"descend")';
    canFire    = ismember(pp.ns,nFire);
    pp.canFire = canFire;
    pp.isConj  = ismember(pp.ns,nConj);
    pp.t       = pp.ns*pp.T;                                                             % [-] (1,N) Maneuver time before TCA
    pp.N       = length(pp.ns);
    n_man      = pp.nMans;
    pp.n_man   = n_man;
    u          = zeros(m,n_man);
    scale      = ones(m,n_man);
    ctrl       = nan(3,n_man);
end
[lim,smdLim,coeffPoC,timeSubtr,xBall,x0,metric] = propDAPoly(pp.DAorder,u,scale,0,0,pp);

%% Optimization
switch pp.solvingMethod
    case 'greedy'
        yF = computeCtrlGreedy(lim,metric,coeffPoC,u, ...
                                   pp.DAorder,scale,n_man);
    case 'global'  % working badly
        yF = computeDvGlobal(lim,coeffPoC,n_man,m,scale);

    case 'nlp'
        yF = computeCtrlNlp(lim,coeffPoC,u,n_man,m,scale);

    case 'moment-relaxation' %not working
        yF = computeDvMomRel(lim,coeffPoC,n_man,m,scale);
end

if pp.fixedDir                                                                  % [-] (3,n_man) Build control matrix node-wise in the case of fixed direction
    ctrl = yF.*pp.thrustDirections(:,1:n_man);
elseif pp.fixedMag                                                              % [-] (3,n_man) Build control matrix node-wise in the case of fixed magnitude
    for j = 1:n_man
        ctrl(1,j) = pp.thrustMagnitude*cos(yF(1,j));
        ctrl(2,j) = pp.thrustMagnitude*sin(yF(1,j))*cos(yF(2,j));
        ctrl(3,j) = pp.thrustMagnitude*sin(yF(1,j))*sin(yF(2,j));
    end
else
    ctrl = yF;                                                                  % [-] (3,n_man) Build control matrix node-wise in the general case
end
simTime = toc - timeSubtr - timeSubtr1;     


%% Validation
metricValPoly = eval_poly(coeffPoC.C,coeffPoC.E,reshape(yF./scale,1,[]),pp.DAorder);

[~,~,~,~,x] = propDAPoly(1,ctrl,scale,1,0,pp);
if pp.metricFlag  == 0 || pp.metricFlag  == 1
    metricValPoly = 10^metricValPoly;
    lim           = 10^lim;
end

%% PostProcess
postProcess(smdLim,xBall,x,lim,ctrl,simTime,pp)
