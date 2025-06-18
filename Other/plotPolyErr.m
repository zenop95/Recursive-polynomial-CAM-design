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
set(groot,'defaultAxesTickLabelInterpreter','latex');  
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
%% User-defined inputs (modifiable)
multiple = 3;                                                                   % [-]     (1,1) flag to activate multiple encounters test case
cislunar = 0;                                                                   % [-]     (1,1) flag to activate cislunar test case
pp = initOpt(multiple,cislunar,1);                                              % [struc] (1,1) Initialize paramters structure with conjunction data
returnTime = -3;                                                                % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
% fireTimes  = [0.6 0.4];                                                               % [-] Example of bi-impulsive maneuvers
fireTimes  = [2.5];                                                               % [-] Example of bi-impulsive maneuvers
% fireTimes  = [0.6, 0.4, -0.6 ,-0.4, -1.6 ,-1.4];                                                               % [-] Example of bi-impulsive maneuvers
% fireTimes = [3.5,2.5,1.5,0.5];                                                  % [-] Example of bi-impulsive maneuvers
% fireTimes = linspace(0.4,0.6,2);                                              % [-] Example of single low-thrust arc
% fireTimes = [linspace(0.4,0.6,2) -linspace(.4,.6,2) -linspace(1.8,2,2)];                        % [-] Example of two low-thrust arcs with different discretization points
pp.cislunar = cislunar;
pp          = defineParams(pp,fireTimes,returnTime);                            % [-] (1,1) Include optimization paramters to parameters structure
% pp.PoCLim   = pp.PoCLim/max(multiple,1);
%% Non-user defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
if pp.flagStability && pp.cislunar
    [CGDir,timeSubtr0] = Cauchy_Green_prop(1,u,pp);
    pp.fixedDir         = true;
    pp.thrustDirections = CGDir;
end
timeSubtr1 = 0;
tic
% If the filtering routine is adpoted, first perform a first-order
% propagation to find the most sensitive maneuvering times
if pp.filterMans
    [~,~,coeff,~,timeSubtr1] = propDA(1,u,scale,0,0,pp);
    gradVec = buildDAArray(coeff.C,coeff.E,1);
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
    ctrl       = nan(3,n_man);
end
aa=tic;
% Propagate the primary orbit and get the PoC coefficient and the position at each TCA

[lim,coeff,timeSubtr,xBall,xRetBall] = propDA(pp.DAorder,u,0,pp);
% if ~pp.flagPoCTot && multiple > 1
%     coeff(pp.n_conj+1) = [];
%     pp.n_constr = pp.n_constr - 1;
% elseif pp.flagPoCTot && multiple > 1
%     coeff(1:pp.n_conj) = [];
%     pp.n_constr = pp.n_constr - pp.n_conj;
% end
metric = coeff(1).C(1);

%% Validation
% nn=201;

Ts = -30:1:30; %[mm/s]
Rs = -100:10:100;
k = 1;
for j = 1:length(Rs)
    j
    for i = 1:length(Ts)
        for k = 1:pp.n_conj
            ctrl = [Rs(j); Ts(i); 0]/1e6/pp.Vsc;
            [~,~,~,x,~,x_sec,deltaTca] = propDA(1,ctrl,1,pp);                      % Validate the solution by forward propagating and computing the real PoC
            x_s    = x_sec;
            STMp   = CWStateTransition(pp.primary.n^(3/2),deltaTca(k)/pp.Tsc,0,1);
            STMs   = CWStateTransition(pp.secondary(k).n^(3/2),deltaTca(k)/pp.Tsc,0,1);
            Cpprop = STMp*pp.Cp(:,:,k)*STMp';
            Csprop = STMs*pp.Cs(:,:,k)*STMs';
            Pp     = Cpprop(1:3,1:3);
            Ps     = Csprop(1:3,1:3);
            r2ep   = rtn2eci(x(1:3,k),x(4:6,k));
            r2es   = rtn2eci(x_s(1:3,k),x_s(4:6,k));
            P      = r2ep*Pp*r2ep' + r2es*Ps*r2es';
            [PB,p,smd] = Bplane(x(:,k),x_s(:,k),P);
            PCs(k) =  poc_Chan(pp.HBR(k),PB,smd,3);
        end
        PoC(j,i) = PoCTot(PCs);                              % [-] (1,1) PoC computed with Chan's formula
    end
end
% load('PoC_comp.mat')
%%
% ts = -50:1:50; %[mm/s]
% rs = -50:1:50;
ts = Ts; %[mm/s]
rs = Rs;

for kk = 1:4
    for j = 1:length(rs)
        for i = 1:length(ts)
            ctrl = [rs(j); ts(i); 0]/1e6/pp.Vsc;
            metricValPoly = eval_poly(coeff.C,coeff.E,ctrl',kk);
            PoCPoly(j,i,kk) = 10^metricValPoly;
        end
    end
end
tt = Ts;%linspace(Ts(1),Ts(end),1e4);
interPoc2 = interp1(Ts,PoC,tt);
figure
for kk = 1:3
    interPoc = interp1(Ts,PoCPoly(:,:,kk),tt);
    plot(tt,abs(interPoc./interPoc2-1))
    hold on
end
legend
hold off

figure
semilogy(tt,abs(interPoc2),'color',[0.5 0.5 0.5],'LineWidth',3)
hold on
for kk = 1:3
    interPoc = interp1(Ts,PoCPoly(:,:,kk),tt);
    semilogy(tt,abs(interPoc))
end
ylim([1e-10,1])
legend
hold off
%%
close all

%%
% for j =1:5
%     err = log10(abs(PoCPoly(:,:,j)));
%     placeFigure
%     s = surfc(X,Y,(err));
%     % clim([-4,0])
%     % view(2)
%     % hold on
%     % s.EdgeColor = 'none';
%     colorbar
% end
% xlabel('R [mm/s]')
% ylabel('T [mm/s]')
%%
[X,Y] = ndgrid(Rs,Ts);
[x,y] = ndgrid(rs,ts);
figure
PoC(PoC>1e-4) = nan;
PoCPoly(PoCPoly>1e-4) = nan;
s = surf(X,Y,log10(PoC));
hold on
s.EdgeColor = 'none';
s1 = surf(x,y,log10(PoCPoly(:,:,2)));
s1.EdgeColor = 'none';
% plot3(zeros(2,1),zeros(2,1),[-20,0],'-o','Linewidth',4)
xlabel('R [mm/s]')
ylabel('T [mm/s]')
colorbar
%%
for j = 1:3
    err = log10(abs(PoCPoly(:,:,j)./PoC-1));
    % err(err<1e-10) = 1e-10;
    err(PoC<1e-12) = nan;
    err(err>0) =nan;
    placeFigure
    % [C,h] = contour(X,Y,err);
    % clabel(C,h)
    s = surfc(X,Y,err);
    clim([-7,0])
    % view(2)
    % hold on
    % s.EdgeColor = 'none';
    colorbar
    xlabel('$\Delta v_R$ [mm/s]')
    ylabel('$\Delta v_T$ [mm/s]')
    zlabel('PoC relative error [-]')
end