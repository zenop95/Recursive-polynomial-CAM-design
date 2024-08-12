beep off
format longG
close all
clear
addpath(genpath('.\data'))
addpath(genpath('.\Functions'))
addpath(genpath('.\CppExec'))
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
%% load variables
n = 2170;
PoCNlp      = nan(6,n);
PoCRec      = nan(6,n);
simTimeNlp  = nan(6,n);
simTimeRec  = nan(6,n);
normNlp       = nan(6,n);
normRec       = nan(6,n);
inPlaneNlp  = nan(6,n);
inPlaneRec  = nan(6,n);
outPlaneNlp = nan(6,n);
outPlaneRec = nan(6,n);
for kk = 2:7
    load(['SimOutput\singleImpulseEsaCases\nlp',num2str(kk),'.mat'])
    % load(['SimOutput\singleImpulseEsaCases\nlp',num2str(kk),'_J2.mat'])
%     load(['SimOutput\doubleImpulseEsaCases\nlp',num2str(kk),'.mat'])
    PoC(PoC>prctile(PoC,99)) = nan;
    PoC(PoC<prctile(PoC,1)) = nan;
    PoCNlp(kk-1,:)     = PoC;
    simTimeNlp(kk-1,:) = simTime;
    dvn = normOfVec(dvs);
    dvn(dvn>prctile(dvn,95)) = nan;
    dvNlp(kk-1,:) = dvn;
    load(['SimOutput\rec',num2str(kk),'.mat'])
    % load(['SimOutput\singleImpulseEsaCases\rec',num2str(kk),'_J2.mat'])
    % load(['SimOutput\doubleImpulseEsaCases\rec',num2str(kk),'.mat'])
    PoC(PoC>prctile(PoC,99)) = nan;
    PoC(PoC<prctile(PoC,1)) = nan;
    PoCRec(kk-1,:)     = PoC;
    simTimeRec(kk-1,:) = simTime;
    dvn = normOfVec(squeeze(dvs));
    dvn(dvn>prctile(dvn,95)) = nan;
    dvRec(kk-1,:) = dvn;
end
% clearvars -except PoCConv simTimeConv PoCNlp simTimeNlp PoCRec simTimeRec dvNlp dvRec n
alsoNlp = 1;
%% Violin plots
figure
A = violin(PoCRec',2:7);
A.LineColor   = [0 0.4470 0.7410];
A.violinColor = [0 0.4470 0.7410];
hold on
if alsoNlp
    B = violin(PoCNlp',2:7);
    B.LineColor = [0.4940 0.1840 0.5560];
    B.violinColor = [0.4940 0.1840 0.5560];
end

hold off
xlabel('Expansion order [-]')
ylabel('PoC [-]')

a = 99.9;
for kk = 2:7
    load(['SimOutput\singleImpulseEsaCases\rec',num2str(kk),'.mat'])
    tcaRec(kk-1,:)     = tcaNewDelta;
    tcaRec(tcaRec>prctile(tcaRec,95)) = nan;
    tcaRec(tcaRec<prctile(tcaRec,5)) = nan;
    inPlaneRec(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
    outPlaneRec(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
    dvn = normOfVec(dvs);
    dvn(dvn>prctile(dvn,a)) = nan;
    normRec(kk-1,:)     = normOfVec(dvs);
    load(['SimOutput\singleImpulseEsaCases\nlp',num2str(kk),'.mat'])
    tcaNlp(kk-1,:)     = tcaNewDelta;
    inPlaneNlp(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
    outPlaneNlp(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
    dvn = normOfVec(dvs);
    dvn(dvn>prctile(dvn,a)) = nan;
    normNlp(kk-1,:)     = normOfVec(dvs);
end
a = 1;

figure
A = violin(tcaRec',2:7);
A.LineColor   = [0 0.4470 0.7410];
A.violinColor = [0 0.4470 0.7410];
% hold on
% if alsoNlp
% B = violin(PoCNlp',2:7);
% end
xlabel('Expansion order [-]')
ylabel('tca diff [s]')

%% Define colors
col1 = [0.9290 0.6940 0.1250];
col2 = [0.4940 0.1840 0.5560];
col3 = [0.4660 0.6740 0.1880];
col4 = [0 0.4470 0.7410];
colorsRec = [linspace(col3(1),col4(1),6)', linspace(col3(2),col4(2),6)', linspace(col3(3),col4(3),6)'];
colorsNlp = [linspace(col1(1),col2(1),6)', linspace(col1(2),col2(2),6)', linspace(col1(3),col2(3),6)'];
colors = [colorsRec;colorsNlp];

%% Delta V histograms
figure()
edgesInP = -8:1.6:8;
edgesOutP = -0.1:.02:0.1;
edgesNorm = 0:20:200;
colororder([colors(1,:);colors(end,:)])
subplot(1,3,1)
histograms([normRec(5,:);normNlp(5,:)],edgesNorm)
xlabel('$\Delta v$ [mm/s]')
xticks([20,100,200])
legend('Recursive','fmincon','interpreter','latex','box','off','Orientation','horizontal')
subplot(1,3,2)
histograms([inPlaneRec(5,:);inPlaneNlp(5,:)],edgesInP)
xlabel('In-plane [deg]')
yticklabels([])
ylim([0,63])
subplot(1,3,3)
histograms([outPlaneRec(5,:);outPlaneNlp(5,:)],edgesOutP)
xlabel('Out-of-plane [deg]')
yticks(0:10:60)
yticklabels([])
ylim([0,63])

%% PoC histograms
figure()
colororder([colors(1,:);colors(end,:)])
edges = (9.9:.0005:10.1)*1e-7;
histograms([PoCRec(5,:);PoCNlp(5,:)],edges)
xlabel('$PoC$ [-]')
