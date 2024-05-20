%% Define path and set figure properties
beep off
format longG
close all
clear
addpath(genpath('.\data'))
addpath(genpath('.\Functions'))
addpath(genpath('.\CppExec'))
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

%% User-defined inputs (modifiable)
multiple = 0;                                                                   % [-]     (1,1) flag to activate multiple encounters test case
cislunar = 0;                                                                   % [-]     (1,1) flag to activate cislunar test case
pp = initOpt(multiple,cislunar,1);                                              % [struc] (1,1) Initialize paramters structure with conjunction data
% fireTimes = [4.4 4.6,3.4 3.6,2.4 2.6,1.4 1.6,0.4 0.6];                                                              % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
fireTimes = 2.5;                                                          % [-] Example of bi-impulsive maneuvers
% fireTimes = [0.49,0.51];                                                          % [-] Example of bi-impulsive maneuvers
% fireTimes = linspace(2.4,2.6,2);                                              % [-] Example of single low-thrust arc
% fireTimes = [linspace(1.4,1.6,3) linspace(2.4,2.6,2)];                        % [-] Example of two low-thrust arcs with different discretization points
pp.cislunar = cislunar;
pp = defineParams(pp,fireTimes);                                                % [-] (1,1) Include optimization paramters to parameters structure

%% Non-user defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N  = pp.N;                                                                      % [-] (1,1) Number of nodes in the propagation
n_man = pp.n_man;                                                               % [-] (1,1) Number of nodes where the maneuver can be performed
if N == 1 && pp.lowThrust; error(['The algorithm needs ' ...
        'at least two nodes to define the low-thrust window']); end
if pp.fixedMag && pp.fixedDir; error(['Both magnitude and direction ' ...
        '                  cannot be fixed at the same time']); end

m     = 3 - pp.fixedMag - 2*pp.fixedDir;  pp.m  = m;                            % [-] (1,1)  Number of optimization variables per node
u     = zeros(m,n_man);                                                         % [-] (m,N)  Ctrl of the unperturbed trajectory
scale = ones(m,n_man);                                                          % [-] (m,N)  ~1 if polynomial scaling is used
ctrl  = nan(3,n_man);                                                           % [-] (3,N)  Initialized ctrl of the optimized trajectory

%% Propagation
timeSubtr1 = 0;
tic
% If the filtering routine is adpoted, first perform a first-order
% propagation to find the most sensitive maneuvering times
if pp.filterMans
    [~,~,coeffPoC,timeSubtr1] = propDA(1,u,scale,0,pp);
    gradVec = buildDAArray(coeffPoC.C,coeffPoC.E,1);
    for j = 1:n_man
        grads(j) = norm(gradVec(1+m*(j-1):m*j));
    end
    [~,thrustNode] = sort(grads,'descend');                                     % Rank the nodes
    thrustNode = thrustNode(1:pp.nMans);                                        % Only keep the first nMans nodes
    % Redefine the problem parameters according to the new nodes definition
    fireTimes      = pp.ns(thrustNode)';
    nConj     = -pp.tca_sep;
    pp.ns      = sort(unique([fireTimes, nConj]),"descend")';
    canFire    = ismember(pp.ns,fireTimes);
    pp.canFire = canFire;
    pp.isConj  = ismember(pp.ns,nConj);
    pp.t       = pp.ns*pp.T;
    pp.N       = length(pp.ns);
    n_man      = pp.nMans;
    pp.n_man   = n_man;
    u          = zeros(m,n_man);
    scale      = ones(m,n_man);
    ctrl       = nan(3,n_man);
end
aa=tic;
% Propagate the primary orbit and get the PoC coefficient and the position at each TCA
[lim,coeffPoC,timeSubtr,xBall,metric] = propDA(pp.DAorder,u,scale,0,pp);
timeProp = toc(aa)-timeSubtr;
%% Optimization
switch pp.solvingMethod
    case 'recursive'
        yF = computeCtrlGreedy(lim,metric,coeffPoC,u, ...
                                   pp.DAorder,scale,n_man);
    case 'recursiveLagrange'
        yF = computeCtrlLagrange(lim,metric,coeffPoC,u, ...
                                   pp.DAorder,scale,n_man,m);
    case 'fmincon'
        yF = computeCtrlNlp(lim,coeffPoC,u,n_man,m,scale);
%     case 'global'  % working badly
%         yF = computeDvGlobal(lim,coeffPoC,n_man,m,scale);
% 
%     case 'moment-relaxation' %not working
%         yF = computeDvMomRel(lim,coeffPoC,n_man,m,scale);
    otherwise 
        error('The solving method should be either recursive or fmincon')
end

if pp.fixedDir                                                                  % [-] (3,n_man) Build control matrix node-wise in the case of fixed direction
    ctrl = yF.*pp.thrustDirections(:,1:n_man);
elseif pp.fixedMag                                                              % [-] (3,n_man) Build control matrix node-wise in the case of fixed magnitude
    for j = 1:n_man
%         ctrl(1,j) = pp.thrustMagnitude*cos(yF(1,j))*sin(yF(2,j));
%         ctrl(2,j) = -pp.thrustMagnitude*cos(yF(1,j))*cos(yF(2,j));
%         ctrl(3,j) = pp.thrustMagnitude*sin(yF(1,j));
        ctrl(1,j) = yF(1,j);
        ctrl(2,j) = yF(2,j);
        ctrl(3,j) = sqrt(pp.thrustMagnitude^2-ctrl(1,j)^2-ctrl(2,j)^2);
    end
else
    ctrl = yF;                                                                  % [-] (3,n_man) Build control matrix node-wise in the general case
end
simTime = toc - timeSubtr - timeSubtr1;     


%% Validation
% metricValPoly = eval_poly(coeffPoC.C,coeffPoC.E,reshape(yF./scale,1,[]), ...    
%                             pp.DAorder);
% metricValPoly = 10^metricValPoly;

[~,~,~,x] = propDA(1,ctrl,scale,1,pp);                                      % Validate the solution by forward propagating and computing the real PoC
lim       = 10^lim;

%% PostProcess
postProcess(xBall,x,lim,ctrl,simTime,pp)
