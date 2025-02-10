beep off
format longG
% close all
% clear
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
name = ["pc4.mat","pc5.mat","pc6.mat","pc7.mat"];
for kk = 1:4
    load("./data/dataConjunctionsESA.mat");
    pclim = [1e-4,1e-5,1e-6,1e-7];
    mdlim = nan(4,1);
    compTime = nan(2170,1);
    xF = nan(6,2170);
    dvs  = nan(3,2170);
    poc  = nan(2170,1);
    md  = nan(2170,1);
    t  = nan(2170,1);
    fracOrb = nan(2170,1);
    for indCase = 1:2170
    t_man = 0.4:0.05:0.9;
    indCase
    multiple = 0;
    returnTime = 0;                                                                 % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
    pp = initOpt(0,0,indCase);
    pp.cislunar = 0;
    pp = defineParams(pp,t_man,returnTime,pclim(kk),mdlim(kk),0);
    pp.nMans   = indCase;                                                            % [bool]   (1,1) Selects how many impulses to use
    N  = pp.N;                                                                      % [-] (1,1) Number of nodes in the propagation
    n_man = pp.n_man;                                                               % [-] (1,1) Number of nodes where the maneuver can be performed
    if N == 1 && pp.lowThrust; error(['The algorithm needs ' ...
            'at least two nodes to define the low-thrust window']); end
    if pp.fixedDir; error(['Both magnitude and direction ' ...
                                'cannot be fixed at the same time']); end
    
    m     = 3 - 2*pp.fixedDir;  pp.m  = m;                            % [-] (1,1)  Number of optimization variables per node
    u     = zeros(m,n_man);                                                         % [-] (m,N)  Ctrl of the unperturbed trajectory
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
    pp.nMans = 1;
    if pp.filterMans
        [~,coeff,timeSubtr1] = propDA(1,u,0,pp);
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
    if ~pp.flagPoCTot && multiple > 1
        coeff(pp.n_conj+1) = [];
        pp.n_constr = pp.n_constr - 1;
    elseif pp.flagPoCTot && multiple > 1
        coeff(1:pp.n_conj) = [];
        pp.n_constr = pp.n_constr - pp.n_conj;
    end
    metric = coeff(1).C(1);
    %% Optimization
    if strcmpi(pp.solvingMethod,'lagrange')
            % [yF,iters,er] = computeCtrlActiveSet(coeff,u,pp);
            yF = computeCtrlRecursive(coeff,u,pp);
    
    elseif strcmpi(pp.solvingMethod,'convex')
            yF = computeCtrlRecursiveConvex(coeff,u,pp);
    
    elseif strcmpi(pp.solvingMethod,'fmincon')
            yF = computeCtrlNlp(coeff,u,pp);
    
    else
        error('The solving method should be either lagrange, fmincon, or convex')
    end
    
    if pp.fixedDir                                                                  % [-] (3,n_man) Build control matrix node-wise in the case of fixed direction
        ctrl = yF.*pp.thrustDirections(:,1:n_man);
    else
        ctrl = yF;                                                                  % [-] (3,n_man) Build control matrix node-wise in the general case
    end
    simTime = toc - timeSubtr - timeSubtr1;     
    ctrl = pp.ctrlMax*ctrl;
    
    %% Validation
    metricValPoly = eval_poly(coeff(1).C,coeff(1).E,reshape(yF,1,[]), ...    
                                pp.DAorder);
    [~,~,~,x,xRetMan,x_sec,deltaTca] = propDA(1,ctrl,1,pp);                      % Validate the solution by forward propagating and computing the real PoC
    if pp.pocType ~= 3
        metricValPoly = 10^metricValPoly;
        lim           = 10^lim;
    end
 
    lim           = 10^lim;
    dvs(:,indCase) = ctrl*pp.Vsc*1e6;
    fracOrb(indCase) = 0.9-0.05*(thrustNode-1);
    t(indCase)     = fracOrb(indCase)*pp.Tsc*pp.T;
    xSec(:,indCase) = x_sec;
    % nodeThrust(:,j)  = thrustNode;
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
    PB(:,:,indCase)     = e2b*P*e2b';
    p      = e2b*(x(1:3)-x_sec(1:3));
    smd    = dot(p,PB(:,:,indCase)\p);
    md(indCase)     = norm(p)*pp.Lsc;
    % PoC(j) = maximumPc(p,PB(:,:,j),pp.HBR);                                        % [-] (1,1) PoC computed with Chan's formula
    poc(indCase) = poc_Chan(pp.HBR,PB(:,:,indCase),smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
    compTime(indCase) = simTime;
    E2B(:,:,indCase) = e2b;
    tcaNewDelta(indCase) = deltaTca;
    xF(:,indCase)   = x.*pp.scaling(1:6);
    end
clearvars -except compTime xF dvs poc md t fracOrb kk name %valErr to dv dto uo retErrV retErrP 
save(name(kk));
end

name = ["md05.mat","md15.mat","md25.mat","md35.mat"];
for kk = 1:4
    load("./data/dataConjunctionsESA.mat");
    mdlim = [0.5,1.5,2.5,3.5];
    pclim = nan(4,1);
    compTime = nan(2170,1);
    xF = nan(6,2170);
    dvs  = nan(3,2170);
    poc  = nan(2170,1);
    md  = nan(2170,1);
    t  = nan(2170,1);
    fracOrb = nan(2170,1);
    for indCase = 1:2170
    t_man = 0.4:0.05:0.9;
    indCase
    multiple = 0;
    returnTime = 0;                                                                 % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
    pp = initOpt(0,0,indCase);
    pp.cislunar = 0;
    pp = defineParams(pp,t_man,returnTime,pclim(kk),mdlim(kk),1);
    pp.nMans   = indCase;                                                            % [bool]   (1,1) Selects how many impulses to use
    N  = pp.N;                                                                      % [-] (1,1) Number of nodes in the propagation
    n_man = pp.n_man;                                                               % [-] (1,1) Number of nodes where the maneuver can be performed
    if N == 1 && pp.lowThrust; error(['The algorithm needs ' ...
            'at least two nodes to define the low-thrust window']); end
    if pp.fixedDir; error(['Both magnitude and direction ' ...
                                'cannot be fixed at the same time']); end
    
    m     = 3 - 2*pp.fixedDir;  pp.m  = m;                            % [-] (1,1)  Number of optimization variables per node
    u     = zeros(m,n_man);                                                         % [-] (m,N)  Ctrl of the unperturbed trajectory
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
    pp.nMans = 1;
    if pp.filterMans
        [~,coeff,timeSubtr1] = propDA(1,u,0,pp);
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
    if ~pp.flagPoCTot && multiple > 1
        coeff(pp.n_conj+1) = [];
        pp.n_constr = pp.n_constr - 1;
    elseif pp.flagPoCTot && multiple > 1
        coeff(1:pp.n_conj) = [];
        pp.n_constr = pp.n_constr - pp.n_conj;
    end
    metric = coeff(1).C(1);
    %% Optimization
    if strcmpi(pp.solvingMethod,'lagrange')
            % [yF,iters,er] = computeCtrlActiveSet(coeff,u,pp);
            yF = computeCtrlRecursive(coeff,u,pp);
    
    elseif strcmpi(pp.solvingMethod,'convex')
            yF = computeCtrlRecursiveConvex(coeff,u,pp);
    
    elseif strcmpi(pp.solvingMethod,'fmincon')
            yF = computeCtrlNlp(coeff,u,pp);
    
    else
        error('The solving method should be either lagrange, fmincon, or convex')
    end
    
    if pp.fixedDir                                                                  % [-] (3,n_man) Build control matrix node-wise in the case of fixed direction
        ctrl = yF.*pp.thrustDirections(:,1:n_man);
    else
        ctrl = yF;                                                                  % [-] (3,n_man) Build control matrix node-wise in the general case
    end
    simTime = toc - timeSubtr - timeSubtr1;     
    ctrl = pp.ctrlMax*ctrl;
    
    %% Validation
    metricValPoly = eval_poly(coeff(1).C,coeff(1).E,reshape(yF,1,[]), ...    
                                pp.DAorder);
    [~,~,~,x,xRetMan,x_sec,deltaTca] = propDA(1,ctrl,1,pp);                      % Validate the solution by forward propagating and computing the real PoC
    if pp.pocType ~= 3
        metricValPoly = 10^metricValPoly;
        lim           = 10^lim;
    end
 
    lim           = 10^lim;
    dvs(:,indCase) = ctrl*pp.Vsc*1e6;
    fracOrb(indCase) = 0.9-0.05*(thrustNode-1);
    t(indCase)     = fracOrb(indCase)*pp.Tsc*pp.T;
    xSec(:,indCase) = x_sec;
    % nodeThrust(:,j)  = thrustNode;
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
    PB(:,:,indCase)     = e2b*P*e2b';
    p      = e2b*(x(1:3)-x_sec(1:3));
    smd    = dot(p,PB(:,:,indCase)\p);
    md(indCase)     = norm(p)*pp.Lsc;
    % PoC(j) = maximumPc(p,PB(:,:,j),pp.HBR);                                        % [-] (1,1) PoC computed with Chan's formula
    poc(indCase) = poc_Chan(pp.HBR,PB(:,:,indCase),smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
    compTime(indCase) = simTime;
    E2B(:,:,indCase) = e2b;
    tcaNewDelta(indCase) = deltaTca;
    xF(:,indCase)   = x.*pp.scaling(1:6);
    end
clearvars -except compTime xF dvs poc md t fracOrb kk name %valErr to dv dto uo retErrV retErrP 
save(name(kk));
end