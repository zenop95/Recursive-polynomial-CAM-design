beep off
format longG
close all
clear
addpath(genpath('.\data'))
addpath(genpath('.\Functions'))
addpath(genpath('.\CppExec'))
addpath(genpath('.\OPM'))
addpath(genpath('.\CDM'))
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
multiple = 0;
for kk = 2:7
for j = 1:2170
pp = initPolyOpt(0,0,j);
pp.cislunar = 0;
% t_man = linspace(0.5,2.5,11);
t_man = 2.5;
pp = definePolyParams(pp,t_man);
pp.DAorder = kk;
N  = pp.N;
n_man = pp.n_man;
if N == 1 && pp.lowThrust; error('The algorithm needs at least two nodes to define the low-thrust window'); end
if pp.fixedMag && pp.fixedDir; error('Both magnitude and direction cannot be fixed at the same time'); end

m     = 3 - pp.fixedMag - 2*pp.fixedDir;   % [-] (1,1)  Number of optimization variables per node
pp.m = m;
u     = zeros(m,n_man);                        % [-] (m,N)  Ctrl of the unperturbed trajectory
scale = ones(m,n_man);                         % [-] (m,N)  ~1 if polynomial scaling is used
dv    = nan(3,n_man);                          % [-] (3,N)  Initialized ctrl of the optimized trajectory

timeSubtr1 = 0;
tic
if pp.filterMans
    [~,~,coeffPoC,timeSubtr1] = propDAPoly(1,u,scale,0,0,pp);
    gradVec = buildDAArray(coeffPoC.C,coeffPoC.E,1);
    for jj = 1:n_man
        grads(jj) = norm(gradVec(1+m*(jj-1):m*jj));
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
switch pp.solvingMethod
    case 'greedy'
        yF = computeCtrlGreedy(lim,metric,coeffPoC,u, ...
                                   pp.DAorder,scale,n_man);
    case 'nlp'
        yF = computeCtrlNlp(lim,coeffPoC,u,n_man,m,scale);
end

if pp.fixedDir
    dv = yF.*pp.thrustDirections(:,1:n_man);
elseif pp.fixedMag
    for j = 1:n_man
        dv(1,j) = pp.thrustMagnitude(j)*cos(yF(1,j));
        dv(2,j) = pp.thrustMagnitude(j)*sin(yF(1,j))*cos(yF(2,j));
        dv(3,j) = pp.thrustMagnitude(j)*sin(yF(1,j))*sin(yF(2,j));
    end
else
    dv = yF;
end
simTime(j) = toc - timeSubtr - timeSubtr1;     


%% Validation
metricValPoly = eval_poly(coeffPoC.C,coeffPoC.E,reshape(yF./scale,1,[]),pp.DAorder);

[~,~,~,~,x] = propDAPoly(1,dv,scale,1,0,pp);
if pp.metricFlag  == 0 || pp.metricFlag  == 1
    metricValPoly = 10^metricValPoly;
    lim           = 10^lim;
end
dvs(:,j) = sum(dv,2)*pp.Vsc*1e6;
xs(:,j)  = x;
% nodeThrust(:,j)  = thrustNode;
e2b    = eci2Bplane(xBall(4:6),pp.x_sTCA(4:6));
e2b    = e2b([1 3],:);
PB     = e2b*pp.P*e2b';
p      = e2b*(x(1:3)-pp.x_sTCA(1:3));
smd    = dot(p,PB\p);
PoC(j) = poc_Chan(pp.HBR,PB,smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
end
save(['simOutput/singleImpulseEsaCasesFixedDir/rec' num2str(pp.DAorder)])
end