beep off
format longG
close all
clear
addpath(genpath('.\data'))
addpath(genpath('.\Functions'))
addpath(genpath('.\CppExec'))
addpath(genpath('.\OPM'))
addpath(genpath('.\CDM'))
addpath(genpath('..\Path'));
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
for j = 1:2170
pp = initPolyOpt(multiple,j);
pp = definePolyParams(pp);
N  = pp.N;
n_man = pp.n_man;
if N == 1 && pp.lowThrust; error('The algorithm needs at least two nodes to define the low-thrust window'); end
if pp.fixedMag && pp.fixedDir; error('Both magnitude and direction cannot be fixed at the same time'); end

m     = 3 - pp.fixedMag - 2*pp.fixedDir;   % [-] (1,1)  Number of optimization variables per node
u     = zeros(m,n_man);                        % [-] (m,N)  Ctrl of the unperturbed trajectory
scale = ones(m,n_man);                         % [-] (m,N)  ~1 if polynomial scaling is used
dv    = nan(3,n_man);                          % [-] (3,N)  Initialized ctrl of the optimized trajectory

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
if pp.singleMan
    [~,~,coeffPoC,timeSubtr1] = propDAPoly(1,u,scale,0,0,pp);
    gradVec  =  buildDAArrayGeneralized(coeffPoC.C,coeffPoC.E,1);
    for j = 1:n_man
        grads(j) = norm(gradVec(1+m*(j-1):m*j));
    end
    [~,thrustNode] = max(grads);
    pp.ns = pp.ns(thrustNode);
    pp.t  = pp.ns*pp.T;                                                              % [-] (1,1) Maneuver time before TCA
    u     = zeros(m,1);
    scale = ones(m,1);
    dv    = nan(3,1);
    pp.N  = 1;
end
[lim,smdLim,coeffPoC,timeSubtr,xBall,x0,metric] = propDAPoly(pp.DAorder,u,scale,0,0,pp);

%% Optimization
switch pp.solvingMethod
    case 'greedy'
        yF = computeDvGreedy(lim,metric,coeffPoC,u, ...
                                   pp.DAorder,pp.maxIter,scale,n_man);
    case 'global' 
        yF = computeDvGlobal(lim,coeffPoC,n_man,m,scale);

    case 'nlp'
        yF = computeDvNlp(lim,coeffPoC,u,n_man,m,scale);

    case 'moment-relaxation' %not working
        yF = computeDvMomRel(lim,coeffPoC,n_man,m,scale);
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
dvs(:,j) = dv;
xs(:,j)  = x;
e2b    = eci2Bplane(xBall(4:6),pp.x_sTCA(4:6));
e2b    = e2b([1 3],:);
PB     = e2b*pp.P*e2b';
p      = e2b*(x(1:3)-pp.x_sTCA(1:3));
smd    = dot(p,PB\p);
PoC(j) = poc_Chan(pp.HBR,PB,smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
end
%%
figure
semilogy(1:10,PoC,'.')
xlabel('Conjunction ID [-]')
ylabel('PoC after maneuver [-]')
grid on

figure
semilogy(1:2170,simTime,'.')
xlabel('Conjunction ID [-]')
ylabel('Computation time [s]')
grid on

% figure
% plot(1:2170,normOfVec(dvs)')
% xlabel('Conjunction ID [-]')
% ylabel('$\Delta v$ [mm/s]')
% grid on

%%
close all
n = 2170;
load('SingleImpulseConvex')
PoC1 = PoC;
simTime1 = simTime;
load('SingleImpulseNlp2')
PoC2 = PoC;
simTime2 = simTime;
load('SingleImpulseNlp5')
PoC3 = PoC;
simTime3 = simTime;
load('SingleImpulseRecursive2.mat')
PoC4 = PoC;
simTime4 = simTime;
load('SingleImpulseRecursive5.mat')
PoC5 = PoC;
simTime5 = simTime;

edges = 0.1:0.1:1.3;
histCount1 = histcounts(simTime1,edges)/n;
histCount2 = histcounts(simTime2,edges)/n;
histCount3 = histcounts(simTime3,edges)/n;
histCount4 = histcounts(simTime4,edges)/n;
histCount5 = histcounts(simTime5,edges)/n;

hists = [histCount1',histCount2',histCount3',histCount4',histCount5']*100;
figure
bar(edges(2:end),hists,'grouped','FaceColor','flat');
grid on
xlabel('Computation time [s]')
legend('Convex','NLP 2nd-order','NLP 5th-order','Recursive 2nd-order','Recursive 5th-order')
ylabel('\% of simulations [-]')


edges = -6.04:0.01:-5.96;
histCount1 = histcounts(log10(PoC1),edges)/n;
histCount2 = histcounts(log10(PoC2),edges)/n;
histCount3 = histcounts(log10(PoC3),edges)/n;
histCount4 = histcounts(log10(PoC4),edges)/n;
histCount5 = histcounts(log10(PoC5),edges)/n;

hists = [histCount1',histCount2',histCount3',histCount4',histCount5']*100;
figure
bar(edges(2:end),hists,'grouped','FaceColor','flat');
grid on
xlabel('log$_{10}$(PoC) [-]')
legend('Convex','NLP 2nd-order','NLP 5th-order','Recursive 2nd-order','Recursive 5th-order')
ylabel('\% of simulations [-]')
