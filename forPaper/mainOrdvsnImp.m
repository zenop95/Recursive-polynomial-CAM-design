% beep off
% format longG
% close all
% clear
% addpath(genpath('.\data'))
% addpath(genpath('.\Functions'))
% addpath(genpath('.\CppExec'))
% addpath(genpath('.\OPM'))
% addpath(genpath('.\CDM'))
% addpath(genpath('.\Path'));
% addpath(genpath('C:\Program Files\Mosek\10.0\toolbox\r2017a'));
% set(0,'DefaultTextInterpreter','latex');
% set(0,'DefaultAxesFontSize',16);
% set(0,'DefaultAxesFontName','Times');
% set(0,'DefaultUicontrolFontName','Times', 'DefaultUicontrolFontSize', 14);
% set(0,'DefaultUitableFontName','Times', 'DefaultUitableFontSize', 14);
% set(0,'DefaultTextFontName','Times', 'DefaultTextFontSize', 14);
% set(0,'DefaultUipanelFontName','Times', 'DefaultUipanelFontSize', 14);
% set(0, 'DefaultLineLineWidth', 1);
% set(0,'defaultfigurecolor',[1 1 1])

%% Initialization variables
multiple = 0;
ts = linspace(0.5,4,8); ts = flip(ts);
for k = 5
for j = 1:8
try
    pp = initPolyOpt(0,0,1060);
    pp.cislunar = 0;
%     t_man = linspace(2.47,2.53,j+1);
    t_man = ts(1:j);
    % t_man = 2.5;
    pp = definePolyParams(pp,t_man);
    pp.DAorder = k+1;
    N  = pp.N;
    n_man = pp.n_man;
    if N == 1 && pp.lowThrust; error('The algorithm needs at least two nodes to define the low-thrust window'); end
    if pp.fixedMag && pp.fixedDir; error('Both magnitude and direction cannot be fixed at the same time'); end
    
    m     = 3 - pp.fixedMag - 2*pp.fixedDir;   % [-] (1,1)  Number of optimization variables per node
    pp.m = m;
    u     = zeros(m,n_man);                        % [-] (m,N)  Ctrl of the unperturbed trajectory
    scale = ones(m,n_man);                         % [-] (m,N)  ~1 if polynomial scaling is used
    ctrl    = nan(3,n_man);                          % [-] (3,N)  Initialized ctrl of the optimized trajectory
    
    timeSubtr1 = 0;
    tic
    if pp.filterMans
        [~,~,coeffPoC,timeSubtr1] = propDA(1,u,scale,0,0,pp);
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
    [lim,smdLim,coeffPoC,timeSubtr,xBall,metric] = propDA(pp.DAorder,u,scale,0,0,pp);
    switch pp.solvingMethod
        case 'greedy'
            yF = computeCtrlGreedy(lim,metric,coeffPoC,u, ...
                                       pp.DAorder,scale,n_man);
        case 'nlp'
            yF = computeCtrlNlp(lim,coeffPoC,u,n_man,m,scale);
    end
    
    if pp.fixedDir
        ctrl = yF.*pp.thrustDirections(:,1:n_man);
    elseif pp.fixedMag
        for j = 1:n_man
            ctrl(1,j) = pp.thrustMagnitude(j)*cos(yF(1,j));
            ctrl(2,j) = pp.thrustMagnitude(j)*sin(yF(1,j))*cos(yF(2,j));
            ctrl(3,j) = pp.thrustMagnitude(j)*sin(yF(1,j))*sin(yF(2,j));
        end
    else
        ctrl = yF;
    end
    simTime(k,j) = toc - timeSubtr - timeSubtr1;     
    
    
    %% Validation
    metricValPoly = eval_poly(coeffPoC.C,coeffPoC.E,reshape(yF./scale,1,[]),pp.DAorder);
    
    [~,~,~,~,x] = propDA(1,ctrl,scale,1,0,pp);
    if pp.metricFlag  == 0 || pp.metricFlag  == 1
        metricValPoly = 10^metricValPoly;
        lim           = 10^lim;
    end
%     t_lt            = pp.t; 
%     t_lt(pp.isConj) = [];
%     dt_lt           = diff(t_lt); 
%     dt_lt           = abs(dt_lt(pp.canFire))';
%     dvs(:,k,j) = sum(ctrl.*dt_lt*pp.Vsc*1e6,2);
    dvs(:,k,j) = sum(ctrl*pp.Vsc*1e6,2);
    xs(:,k,j)  = x;
    % nodeThrust(:,j)  = thrustNode;
    e2b    = eci2Bplane(xBall(4:6),pp.x_sTCA(4:6));
    e2b    = e2b([1 3],:);
    PB     = e2b*pp.P*e2b';
    p      = e2b*(x(1:3)-pp.x_sTCA(1:3));
    smd    = dot(p,PB\p);
    PoC(k,j) = poc_Chan(pp.HBR,PB,smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
    % nodeThrust(:,j) = thrustNode;
catch
    simTime(k,j) = nan;
    dvs(:,k,j) = nan(3,1);
    xs(:,k,j)  = nan(6,1);
    PoC(k,j) = nan;
end
end
end
% save('simOutput/ordVsNImp/recFixedDir')

%%
kk = 5;
jj = 8;
for j = 1:jj
    for k = 1:kk
        dvNorm(k,j) = normOfVec(dvs(:,k,j));
    end
end
PoC(PoC == 0) = nan;
dvNorm(dvNorm == 0) = nan;
simTime(simTime == 0) = nan;
[x,y] = meshgrid(1:jj,2:kk+1);

figure
subplot(3,1,1)
% surf(x,y,PoC(1:kk,1:jj))
s = scatter(reshape(x,1,[]),reshape(y,1,[]),[],reshape(abs(PoC(1:kk,1:jj)-1e-6),1,[]),'filled');
s.SizeData = 100;
hold on
hold off
% xlabel('Nodes number [-]')
xticklabels([])
ylabel('$n$ [-]')
grid on
box on
% set(gca, 'YDir','reverse')
% hXLabel = get(gca,'XLabel');
% set(hXLabel,'rotation',75,'VerticalAlignment','middle')
cb = colorbar;
ylabel(cb,'PoC Error [-]','Interpreter','latex')
% caxis([8,10]*1e-7)
view(2)
xticks(1:jj)
yticks(2:kk+1)
xlim([0,jj+1])
ylim([1,kk+2])

subplot(3,1,2)
% surf(x,y,dvNorm(1:kk,1:jj))
s = scatter(reshape(x,1,[]),reshape(y,1,[]),[],reshape(dvNorm(1:kk,1:jj),1,[]),'filled');
s.SizeData = 100;
ylabel('$n$ [-]')
grid on
box on
xticklabels([])
% set(gca, 'YDir','reverse')
% hXLabel = get(gca,'XLabel');
% set(hXLabel,'rotation',75,'VerticalAlignment','middle')
cb = colorbar;
ylabel(cb,'$||\Delta v||$ [mm/s]','Interpreter','latex')
view(2)
xticks(1:jj)
yticks(2:kk+1)
xlim([0,jj+1])
ylim([1,kk+2])

subplot(3,1,3)
% surf(x,y,log10(simTime(1:kk,1:jj)))
s = scatter(reshape(x,1,[]),reshape(y,1,[]),[],reshape(log10(simTime(1:kk,1:jj)),1,[]),'filled');
s.SizeData = 100;
% xlabel('Nodes number [-]')
ylabel('$n$ [-]')
grid on
box on
% set(gca, 'YDir','reverse')
% hXLabel = get(gca,'XLabel');
% set(hXLabel,'rotation',75,'VerticalAlignment','middle')
cb = colorbar;
view(2)
xlabel('Number of thrusting nodes [-]')
ylabel(cb,'Time [log$_{10}$(s)]','Interpreter','latex')
xticks(1:jj)
yticks(2:kk+1)
xlim([0,jj+1])
ylim([1,kk+2])
