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
t_man = [0.5, -0.5 0];
returnTime = -1;                                                           % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
for kk = 5
for j = 1:2170
j
pp = initOpt(0,0,j);
pp.cislunar = 0;
pp = defineParams(pp,t_man,returnTime);
pp.DAorder = kk;
pp.nMans   = j;                                                            % [bool]   (1,1) Selects how many impulses to use
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
[lim,coeff,timeSubtr,xBall] = propDA(pp.DAorder,u,scale,0,pp);
if ~pp.flagPoCTot && multiple > 1
    coeff(pp.n_conj+1) = [];
    pp.n_constr = pp.n_constr - 1;
end
metric = coeff(1).C(1);
%% Optimization
if any(strcmpi(pp.solvingMethod,{'lagrange','convex','newton'}))
        yF = computeCtrlRecursive(coeff,u,scale,pp);

elseif strcmpi(pp.solvingMethod,'fmincon')
        yF = computeCtrlNlp(coeff,u,scale,pp);

else
    error('The solving method should be either lagrange, fmincon, or convex')
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
% metricValPoly = eval_poly(coeff.C,coeff.E,reshape(yF./scale,1,[]),pp.DAorder);

[~,~,~,x,xRet0] = propDA(1,ctrl,scale,1,pp);
errRetEci = xRet0 - pp.xReference;
errP(j) = norm(errRetEci(1:3))*pp.Lsc*1e3;
errV(j) = norm(errRetEci(4:6))*pp.Vsc*1e6;
% metricValPoly = 10^metricValPoly;
lim           = 10^lim;
dvs(:,:,j) = ctrl*pp.Vsc*1e6;
xs(:,j)  = x;
xSec(:,j) = pp.x_sTCA;
% nodeThrust(:,j)  = thrustNode;
e2b    = eci2Bplane(xBall(4:6),pp.x_sTCA(4:6));
e2b    = e2b([1 3],:);
PB(:,:,j)     = e2b*pp.P*e2b';
p      = e2b*(x(1:3)-pp.x_sTCA(1:3));
smd    = dot(p,PB\p);
PoC(j) = poc_Chan(pp.HBR,PB(:,:,j),smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
compTime(j) = simTime;
E2B(:,:,j) = e2b;
% nodeThrust(:,j) = thrustNode;
end
clearvars -except errP errV dvs xs PoC compTime PB E2B xSec pp
save(['simOutput/IAC/rec' num2str(pp.DAorder)])
end