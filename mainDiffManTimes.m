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
tMan     = linspace(5.5,0.05,nMan);
for kk = 2:4
for j = 1:nMan
pp = initOpt(0,0,1);
pp.cislunar = 0;
pp = defineParams(pp,tMan(j),0);
pp.cislunar = 0;
pp.nMans   = 1;                                                            % [bool]   (1,1) Selects how many impulses to use
j
pp.DAorder = kk;
N  = pp.N;                                                                      % [-] (1,1) Number of nodes in the propagation
n_man = pp.n_man;                                                               % [-] (1,1) Number of nodes where the maneuver can be performed
if N == 1 && pp.lowThrust; error(['The algorithm needs ' ...
        'at least two nodes to define the low-thrust window']); end
if pp.fixedMag && pp.fixedDir; error(['Both magnitude and direction ' ...
                            'cannot be fixed at the same time']); end

m     = 3 - pp.fixedMag - 2*pp.fixedDir;  pp.m  = m;                            % [-] (1,1)  Number of optimization variables per node
u     = zeros(m,n_man);                                                         % [-] (m,N)  Ctrl of the unperturbed trajectory
scale = ones(m,n_man);                                                          % [-] (m,N)  ~1 if polynomial scaling is used
ctrl  = nan(3,n_man);                                                           % [-] (3,N)  Initialized ctrl of the optimized trajectory

%% Propagation
timeSubtr0 = 0;
timeSubtr1 = 0;
tic;
% Propagate the primary orbit and get the PoC coefficient and the position at each TCA

[lim,coeff,timeSubtr,xBall] = propDA(pp.DAorder,u,scale,0,pp);
if ~pp.flagPoCTot && multiple > 1
    coeff(pp.n_conj+1) = [];
    pp.n_constr = pp.n_constr - 1;
elseif pp.flagPoCTot && multiple > 1
    coeff(1:pp.n_conj) = [];
    pp.n_constr = pp.n_constr - pp.n_conj;
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
ctrl = pp.ctrlMax*ctrl;

%% Validation
metricValPoly = eval_poly(coeff(1).C,coeff(1).E,reshape(yF./scale,1,[]), ...    
                            pp.DAorder);
% distValPoly = eval_poly(coeff(2).C,coeff(2).E,reshape(yF./scale,1,[]), ...    
                            % pp.DAorder)*pp.Lsc;
[~,~,~,x,xRet0,x_sec,deltaTca] = propDA(1,ctrl,scale,1,pp);                      % Validate the solution by forward propagating and computing the real PoC
if pp.pocType ~= 3
    metricValPoly = 10^metricValPoly;
    lim           = 10^lim;
end
lim           = 10^lim;
dvs(:,:,j) = squeeze(ctrl*pp.Vsc*1e6);
xs(:,j)  = x;
xSec(:,j) = x_sec;
STMp   = CWStateTransition(pp.primary.n^(3/2),deltaTca/pp.Tsc,0,1);
STMs   = CWStateTransition(pp.secondary.n^(3/2),deltaTca/pp.Tsc,0,1);
Cpprop = STMp*pp.Cp*STMp';
Csprop = STMs*pp.Cs*STMs';
Pp     = Cpprop(1:3,1:3);
Ps     = Csprop(1:3,1:3);
r2ep   = rtn2eci(x(1:3),x(4:6));
r2es   = rtn2eci(x_sec(1:3),x_sec(4:6));
P      = r2ep*Pp*r2ep' + r2es*Ps*r2es';
e2b    = eci2Bplane(x(4:6),x_sec(4:6));
e2b    = e2b([1 3],:);
PB(:,:,j)     = e2b*P*e2b';
p      = e2b*(x(1:3)-x_sec(1:3));
smd    = dot(p,PB(:,:,j)\p);
PoC(j) = poc_Chan(pp.HBR,PB(:,:,j),smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
compTime(j) = simTime;
E2B(:,:,j) = e2b;
tcaNewDelta(j) = deltaTca;                                      % [-] (1,1) PoC computed with Chan's formula
end
clearvars -except dvs xs xSec PoC compTime tcaNewDelta pp t_man xBall tMan
save(['simOutput/diffMan' num2str(pp.DAorder)])
end
%%
nMan = 500;
for jj = 2:5
    load(['SimOutput\diffMan',num2str(jj),'.mat'])
    tca(jj-1,:)     = tcaNewDelta;
    PoCs(jj-1,:)     = PoC;
    simTimes(jj-1,:) = compTime;
    movmeanTime(jj-1,:) = movmean(compTime,100);
    dv = (squeeze(dvs));
end

%%
col1 = [0.56,0.80,0.90];
col2 = [0.00,0.45,0.74];
colors = [linspace(col1(1),col2(1),4)', linspace(col1(2),col2(2),4)', linspace(col1(3),col2(3),4)'];
figure
colororder(colors)
subplot(2,1,1)
scatter(tMan,abs(1e-6-PoCs),5,'filled')
xticklabels([])
ylabel('PoC error [-]')
set(gca, 'XDir','reverse')
grid on
axis tight
legend('$n=2$','$n=3$','$n=4$','$n=5$','interpreter','latex','Orientation','horizontal')
box on

subplot(2,1,2)
scatter(tMan,abs(tca(4,:))',5,'filled')
xlabel('Orbits to TCA [-]')
ylabel('TCA shift [s]')
set(gca, 'XDir','reverse')
grid on
axis tight
box on

figure
plot(tMan,simTimes','.')
xlabel('Orbits to TCA [-]')
ylabel('Computation time [s]')
set(gca, 'XDir','reverse')
grid on

figure
plot(tMan,movmeanTime')
xlabel('Orbits to TCA [-]')
ylabel('Computation time [s]')
set(gca, 'XDir','reverse')
grid on

figure
scatter(tMan,squeeze(dvs)',5,'filled')
hold on
scatter(tMan,normOfVec(squeeze(dvs)),5,'filled')
legend('R','T','N','$||\cdot||$','Interpreter','Latex','Orientation','horizontal')
xlabel('Orbits to TCA [-]')
ylabel('$\Delta v$ [mm/s]')
grid on
set(gca, 'XDir','reverse')
hold off
box on
axis tight

%% Ellipse B-plane
e2b = eci2Bplane(xBall(4:6,1),pp.x_sTCA(4:6));
e2b = e2b([1 3],:);
PB  = e2b*pp.P*e2b';
smdLim   = PoC2SMD(PB, pp.HBR, pp.PoCLim, 3, 1, 1e-3, 200);   % [-] (1,1) SMD limit computed with Chan's formula
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
    pNewB(:,j) = e2b*(xs(1:3,j)-xSec(1:3,j))*pp.Lsc;
end
plot(ellB(2,:),ellB(1,:),'k');
plot(pOldB(2),pOldB(1),'k','marker','diamond','HandleVisibility','off')
s = scatter(pNewB(2,275:321),pNewB(1,275:321),[],tMan(275:321)','filled');
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
xlim([-.4,0])
ylim([0,.15])


simTimes(1,1) = nan;
figure
A = violin(simTimes',2:5);
hold off
xlabel('Expansion order [-]')
