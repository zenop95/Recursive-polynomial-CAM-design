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
nMan     = 500;
tMan     = linspace(4,0.05,nMan);
for j = 1:nMan
pp = initPolyOpt(0,1,1);
pp.cislunar = 1;
pp = definePolyParams(pp,tMan(j));
N  = pp.N;
n_man = pp.n_man;
if N == 1 && pp.lowThrust; error('The algorithm needs at least two nodes to define the low-thrust window'); end
if pp.fixedMag && pp.fixedDir; error('Both magnitude and direction cannot be fixed at the same time'); end

m     = 3 - pp.fixedMag - 2*pp.fixedDir;   % [-] (1,1)  Number of optimization variables per node
pp.m = m;
u     = zeros(m,n_man);                        % [-] (m,N)  Ctrl of the unperturbed trajectory
scale = ones(m,n_man);                         % [-] (m,N)  ~1 if polynomial scaling is used
dv    = nan(3,n_man);                          % [-] (3,N)  Initialized ctrl of the optimized trajectory

%% Propagation
timeSubtr1 = 0;
tic
if pp.filterMans
    [~,~,coeffPoC,timeSubtr1] = propDA(1,u,scale,0,0,pp);
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
[lim,smdLim,coeffPoC,timeSubtr,xBall,metric] = propDA(pp.DAorder,u,scale,0,0,pp);

%% Optimization
switch pp.solvingMethod
    case 'greedy'
        yF = computeCtrlGreedy(lim,metric,coeffPoC,u, ...
                                   pp.DAorder,scale,n_man);
    case 'global' 
        yF = computeDvGlobal(lim,coeffPoC,n_man,m,scale);

    case 'nlp'
        yF = computeCtrlNlp(lim,coeffPoC,u,n_man,m,scale);

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

[~,~,~,~,x] = propDA(1,dv,scale,1,0,pp);
metricValPoly = 10^metricValPoly;
lim           = 10^lim;
dvs(:,j) = dv*pp.Vsc*1e6;
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
semilogy(tMan,PoC)
xlabel('Orbits to TCA [-]')
ylabel('PoC after maneuver [-]')
grid on

figure
semilogy(tMan,simTime,'.')
xlabel('Orbits to TCA [-]')
ylabel('Computation time [s]')
grid on

figure
plot(tMan,dvs')
hold on
plot(tMan,normOfVec(dvs),'k')
legend('R','T','N','$||\cdot||$','Interpreter','Latex')
xlabel('Orbits to TCA [-]')
ylabel('$\Delta v$ [mm/s]')
grid on
set(gca, 'XDir','reverse')

%% Ellipse B-plane
e2b = eci2Bplane(xBall(4:6,1),pp.x_sTCA(4:6));
e2b = e2b([1 3],:);
PB  = e2b*pp.P*e2b';
[semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
a          = semiaxes(1)*pp.Lsc;
b          = semiaxes(2)*pp.Lsc;
tt         = 0:0.001:2*pi;
xx         = a*cos(tt);
yy         = b*sin(tt);
ellCov     = [xx; yy];
ellB       = nan(2,length(tt));
for k = 1:length(tt)
    ellB(:,k) = cov2b*ellCov(:,k);
end
figure
hold on    
pOldB = e2b*(xBall(1:3)-pp.x_sTCA(1:3))*pp.Lsc;
for j = 1:nMan
    pNewB(:,j) = e2b*(xs(1:3,j)-pp.x_sTCA(1:3))*pp.Lsc;
end
plot(ellB(2,:),ellB(1,:),'k');
plot(pOldB(2),pOldB(1),'k','marker','diamond','HandleVisibility','off')
s = scatter(pNewB(2,:),pNewB(1,:),[],tMan','filled');
s.SizeData = 20;
colormap(flipud(jet));
grid on 
xlabel('$\zeta$ [km]')
ylabel('$\xi$ [km]')
hold off
cb = colorbar;
if pp.cislunar; str = 'Days'; else; str = 'Orbits'; end
ylabel(cb, [str, ' to TCA [-]'])
clearvars -except PoC dvs xs pp tMan nMan simTime xBall smdLim
save(['simOutput/diffManTimesordCislunar' num2str(pp.DAorder)])