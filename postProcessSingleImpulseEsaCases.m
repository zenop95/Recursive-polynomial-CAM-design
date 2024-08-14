beep off
format longG
% close all
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
simFolder = 'SimOutput\impulsive\';
n = 2170;
PoCNlp      = nan(6,n);
PoCRec      = nan(6,n);
simTimeNlp  = nan(6,n);
simTimeRec  = nan(6,n);
normNlp     = nan(6,n);
normRec     = nan(6,n);
inPlaneNlp  = nan(6,n);
inPlaneRec  = nan(6,n);
outPlaneNlp = nan(6,n);
outPlaneRec = nan(6,n);
for kk = 2:7
    load([simFolder,'nlp',num2str(kk),'.mat'])
    PoC(PoC>prctile(PoC,95)) = nan;
    PoC(PoC<prctile(PoC,5)) = nan;
    tcaNewDelta(tcaNewDelta>prctile(tcaNewDelta,99.95)) = nan;
    tcaNewDelta(tcaNewDelta<prctile(tcaNewDelta,0.05)) = nan;
    compTime(compTime>prctile(compTime,99)) = nan;
    PoCNlp(kk-1,:)     = PoC;
    tcaRec(kk-1,:)     = tcaNewDelta;
    simTimeNlp(kk-1,:) = compTime;
    inPlaneNlp(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
    outPlaneNlp(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
    dvn = normOfVec(squeeze(dvs));
    dvn(dvn>prctile(dvn,99.9)) = nan;
    normNlp(kk-1,:)     = dvn;
    load([simFolder,'rec',num2str(kk),'.mat'])
    PoC(PoC>prctile(PoC,95)) = nan;
    PoC(PoC<prctile(PoC,5)) = nan;
    tcaNewDelta(tcaNewDelta>prctile(tcaNewDelta,99.95)) = nan;
    tcaNewDelta(tcaNewDelta<prctile(tcaNewDelta,0.05)) = nan;
    tcaNlp(kk-1,:)     = tcaNewDelta;
    compTime(compTime>prctile(compTime,99)) = nan;
    PoCRec(kk-1,:)     = PoC;
    simTimeRec(kk-1,:) = compTime;
    inPlaneRec(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
    outPlaneRec(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
    dvn = normOfVec(squeeze(dvs));
    dvn(dvn>prctile(dvn,99.9)) = nan;
    normRec(kk-1,:)     = dvn;  
end
alsoNlp = 1;
%% Violin plots
figure
A = violin(PoCRec'*1e6-1,2:7);
hold on
% if alsoNlp
%     B = violin(PoCNlp'-1e-6,2:7);
%     B.LineColor = [0.4940 0.1840 0.5560];
%     B.violinColor = [0.4940 0.1840 0.5560];
%     B.mediancolor = 'r';
%     B.meancolor   = 'r';
% end
hold off
xlabel('Expansion order [-]')
ylabel('PoC [-]')

% figure
% subplot(3,1,3)
% A = violin(simTimeRec',2:7);
% hold on
% if alsoNlp
%     B = violin(simTimeNlp',2:7);
%     B.LineColor = [0.4940 0.1840 0.5560];
%     B.violinColor = [0.4940 0.1840 0.5560];
%     B.mediancolor = 'r';
%     B.meancolor   = 'r';
% end
% hold off
% xlabel('Expansion order [-]')
% ylabel('Computation Time [s]')

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

% TCA histogram
figure('Renderer', 'painters', 'Position', [300 300 560 250])
colororder([colors(1,:);colors(end,:)])
edges = -1:0.07:1;
histograms([tcaRec(5,:);tcaNlp(5,:)],edges)
xlabel('$\Delta t_{CA}$ [s]')
legend('Recursive','fmincon','interpreter','latex','box','off','Orientation','horizontal')
