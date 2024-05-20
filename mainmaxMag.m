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
pp = initOpt(multiple,cislunar,1219);                                              % [struc] (1,1) Initialize paramters structure with conjunction data
fireTimes = [4.4 4.6,3.4 3.6,2.4 2.6,1.4 1.6,0.4 0.6];                                                              % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
% fireTimes = 2.5;                                                          % [-] Example of bi-impulsive maneuvers
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
ctrl  = zeros(3,n_man);                                                           % [-] (3,N)  Initialized ctrl of the optimized trajectory

%% Propagation
timeSubtr1 = 0;
tic
% If the filtering routine is adpoted, first perform a first-order
% propagation to find the most sensitive maneuvering times
if pp.filterMans
    [~,coeffPoC,timeSubtr1] = propDA(1,u,scale,0,pp);
    gradVec = buildDAArray(coeffPoC.C,coeffPoC.E,1);
    for j = 1:n_man
        grads(j) = norm(gradVec(1+m*(j-1):m*j));
    end
    [~,thrustNode] = sort(grads,'descend');                                     % Rank the nodes
    thrustNode = thrustNode(1:pp.nMans);                                        % Only keep the first nMans nodes
    % Redefine the problem parameters according to the new nodes definition
    nConj      = -pp.tca_sep;
    if pp.lowThrust
        fireTimes = pp.ns(1:end-1)';
        pp.ns     = sort(unique([fireTimes, nConj]),"descend")';
        canFire   = ismember(pp.ns,fireTimes(1:2:end));
    else
        fireTimes  = pp.ns(thrustNode)';
        pp.ns       = sort(unique([fireTimes, nConj]),"descend")';
        canFire     = ismember(pp.ns,fireTimes);
    end
    pp.canFire  = canFire;
    pp.isConj   = ismember(pp.ns,nConj);
    pp.t        = pp.ns*pp.T;
    pp.N        = length(pp.ns);
    n_man       = pp.nMans;
    pp.n_man    = n_man;
    u           = zeros(m,n_man);
    scale       = ones(m,n_man);
    ctrl        = zeros(3,0);
    pp.canFire0 = pp.canFire;
    pp.isConj0  = pp.isConj;
    pp.t0       = pp.t;
    pp.N0       = pp.N;
    pp.n_man0   = pp.n_man;
end
aa=tic;
it = 0;
e = 1;
pp.t = [pp.t(1:2); pp.t(end)];
scale0 = scale;
% Propagate the primary orbit and get the PoC coefficient and the position at each TCA
while e > 1e-3 && it <= pp.n_man0
pp.t       = [pp.t0(1:it*2+2); pp.t0(end)];
pp.canFire = [pp.canFire0(1:it*2+2); pp.canFire0(end)];
pp.isConj  = [pp.isConj0(1:it*2+2); pp.isConj0(end)];
u = [ctrl(:,1:it), zeros(3,it-size(ctrl(:,1:it),2)+1)];
it = it + 1;
scale = scale0(:,1:it);
pp.N = length(pp.t);
pp.n_man = it;
n_man = pp.n_man;
[lim,coeffPoC,timeSubtr,xBall,metric] = propDA(pp.DAorder,u,scale,0,pp);
c = boolean(ones(length(coeffPoC.C),1));
C = coeffPoC.C; E = coeffPoC.E;
if it > 1
for i = 1:length(coeffPoC.C)
    if any(coeffPoC.E(i,1:end-3))
       c(i) = false;
    end
end
coeffPoC.C = coeffPoC.C(c);
coeffPoC.E = coeffPoC.E(c,end-2:end);
end
timeProp = toc(aa)-timeSubtr;
%% Optimization
switch pp.solvingMethod
    case 'recursive'
        yF = computeCtrlGreedy(lim,metric,coeffPoC,u(:,end), ...
                                   pp.DAorder,scale(:,end),1);
    case 'recursiveLagrange'
        yF = computeCtrlLagrange(lim,metric,coeffPoC,u, ...
                                   pp.DAorder,scale,n_man,m);
    case 'fmincon'
        yF = computeCtrlNlp(lim,coeffPoC,u,n_man,m,scale);
    otherwise 
        error('The solving method should be either recursive or fmincon')
end


ctrl = [ctrl yF];                                                                  % [-] (3,n_man) Build control matrix node-wise in the general case
for j = 1:n_man
    if norm(ctrl(:,j)) > pp.maxMag
        ctrl(:,j) = pp.maxMag*normalize(ctrl(:,j),'norm');
    end
end
simTime = toc - timeSubtr - timeSubtr1;     

[~,~,~,x] = propDA(1,ctrl,scale,1,pp);                                      % Validate the solution by forward propagating and computing the real PoC
xb     = xBall;
x_s    = pp.x_sTCA;
e2b    = eci2Bplane(xb(4:6),x_s(4:6));
e2b    = e2b([1 3],:);
PB     = e2b*pp.P*e2b';
p      = e2b*(x(1:3)-x_s(1:3));
smd    = dot(p,PB\p);
PoC = poc_Chan(pp.HBR,PB,smd,3);                                        % [-] (1,1) PoC computed with Chan's formula

e = log10(PoC) - lim;

end
%% Validation
lim       = 10^lim;

%% PostProcess
postProcess(xBall,x,lim,ctrl,simTime,pp)
