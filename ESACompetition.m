beep off
format longG
clear
addpath(genpath('.\data'))
addpath(genpath('.\Functions'))
addpath(genpath('.\CppExec'))
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

%% Run Simulations
n = 3;
for indCase = 1:n
    disp(['Case ', num2str(indCase)])
    pp.orbit = 'leo';
    pp       = simProperties(pp.orbit,indCase);
    N        = pp.N;
    pp.timeSubtr = 0;
    pp       = defIndConv(pp);
    pp.nUpd  = 0;
    pp       = findFiringNodes(pp);
    pp       = generateAida(pp);    
    tic
    iterSmd              = 0;
    newTraj              = nan(6,N);
    [newTraj,iterSmd,pp] = majorSolve(iterSmd,newTraj,pp);
    pp.simTime = toc - pp.timeSubtr;
    pp = validate(pp);
    N          = pp.N;
    NCA        = pp.NCA;
    simT       = pp.simTime;
    u          = pp.u;
    statesLin  = pp.majorIter(end).x;
    statesVal  = pp.validationAbsTraj;
    try tErr   = norm(pp.validationAbsTraj - pp.cartTarget); catch; end
    pocConst   = pp.PoCConst;
    Dv         = pp.DvTot;
    dvAll      = pp.dv;
    scaling    = pp.scaling;
    t          = pp.t*pp.Tsc/3600; %time in hours
    smdLim     = pp.majorIter(end).sqrMahaLim;
    pOld       = pp.ballisticRelTraj(1:3,NCA);
    pNew       = pp.validationTraj(1:3,NCA);    
    e2b        = pp.e2b([1 3],:);
    PB         = e2b*pp.majorIter(end).P(:,:,NCA)*e2b';
    majorN     = length(pp.majorIter);
    iterErrFin = pp.majorIter(end).err;
    clear pp;
    save(['Base_',num2str(indCase)])
end
%% PostProcess
for k = 1:n
    load(['Base_',num2str(k)]);
    poc(k)        = pocConst;
    dv(k)         = Dv;
    xF(:,k)       = pNew*scaling(1);
    simTime(k)    = simT;
    valErr(k)     = max(normOfVec(statesLin-statesVal));
%         targErr(k)    = tErr;
    nMaj(k)       = majorN;
end
cases = 1:n;

figure
plot(cases,dv)
grid on
xlabel('Case id')
ylabel('$\Delta v$ [m/s]')

figure
plot(cases,simTime)
grid on
xlabel('Case id')
ylabel('$T_{sim}$ [s]')

figure
plot(cases,valErr)
grid on
xlabel('Case id')
ylabel('Validation error [-]')

figure
plot(cases,nMaj)
grid on
xlabel('Case id')
ylabel('Number of Major Iterations [-]')

figure
plot(cases,poc)
grid on
xlabel('Case id')
ylabel('$P_C$ [-]')