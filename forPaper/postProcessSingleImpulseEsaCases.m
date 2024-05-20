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
load('SimOutput\singleImpulseEsaCases\convex.mat')
PoCConv     = PoC;
simTimeConv = simTime;
% for kk = 2:7
% %     load(['SimOutput\singleImpulseEsaCases\nlp',num2str(kk),'.mat'])
%     load(['SimOutput\singleImpulseEsaCases\nlp',num2str(kk),'_J2.mat'])
% %     load(['SimOutput\doubleImpulseEsaCases\nlp',num2str(kk),'.mat'])
%     PoC(PoC>prctile(PoC,99)) = nan;
%     PoC(PoC<prctile(PoC,1)) = nan;
%     PoCNlp(kk-1,:)     = PoC;
%     simTimeNlp(kk-1,:) = simTime;
%     dvn = normOfVec(dvs);
%     dvn(dvn>prctile(dvn,95)) = nan;
%     dvNlp(kk-1,:) = dvn;
% %     load(['SimOutput\singleImpulseEsaCases\rec',num2str(kk),'.mat'])
%     load(['SimOutput\singleImpulseEsaCases\rec',num2str(kk),'_J2.mat'])
% %     load(['SimOutput\doubleImpulseEsaCases\rec',num2str(kk),'.mat'])
%     PoC(PoC>prctile(PoC,99)) = nan;
%     PoC(PoC<prctile(PoC,1)) = nan;
%     PoCRec(kk-1,:)     = PoC;
%     simTimeRec(kk-1,:) = simTime;
%     dvn = normOfVec(dvs);
%     dvn(dvn>prctile(dvn,95)) = nan;
%     dvRec(kk-1,:) = dvn;
% end
% clearvars -except PoCConv simTimeConv PoCNlp simTimeNlp PoCRec simTimeRec dvNlp dvRec n
alsoNlp = 1;
%% Violin plots
% figure
% A = iosr.statistics.boxPlot(2:7,PoCRec');
% A.mediancolor = 'b';
% A.outliersize = nan;
% A.showMean = true;
% A.percentile = [50 50]; 
% A.lineStyle = 'none';
% A.meancolor = 'b';
% A.showViolin = true;
% A.showScatter = false;
% A.LineColor   = [0 0.4470 0.7410];
% A.violinColor = [0 0.4470 0.7410];
% A.violinAlpha = 0.4;
% if alsoNlp
%     hold on
%     B = iosr.statistics.boxPlot(2:7,PoCNlp');
%     B.mediancolor = 'r';
%     B.outliersize = nan;
%     B.showMean = true;
%     B.percentile = [50 50]; 
%     B.lineStyle = 'none';
%     B.meancolor = 'r';
%     B.showViolin = true;
%     B.showScatter = false;
%     B.LineColor = [0.4940 0.1840 0.5560];
%     B.violinColor = [0.4940 0.1840 0.5560];
%     B.violinAlpha = 0.3;
%  end
% xlabel('Expansion order [-]')
% ylabel('PoC [-]')
% grid on
% box on
a = 99.9;
for kk = 2:7
    load(['SimOutput\singleImpulseEsaCases\rec',num2str(kk),'.mat'])
    inPlaneRec(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
    outPlaneRec(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
    dvn = normOfVec(dvs);
    dvn(dvn>prctile(dvn,a)) = nan;
    normRec(kk-1,:)     = normOfVec(dvs);
    load(['SimOutput\singleImpulseEsaCases\nlp',num2str(kk),'.mat'])
    inPlaneNlp(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
    outPlaneNlp(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
    dvn = normOfVec(dvs);
    dvn(dvn>prctile(dvn,a)) = nan;
    normNlp(kk-1,:)     = normOfVec(dvs);
end
a = 1;
% figure
% A = iosr.statistics.boxPlot(inPlaneRec);
% A.mediancolor = 'b';
% A.outliersize = 1;
% A.showMean = true;
% A.percentile = [25 75]; 
% % A.lineStyle = 'none';
% A.meancolor = 'b';
% % A.showViolin = true;
% A.showScatter = true;
% A.scatterLayer = 'bottom';         % The layer of scatter plots with respect to the boxes.
% A.scatterMarker = 'x';
% A.scatterColor = 'g';
% A.scatterSize   = 5;
% % A.showViolin = true;
% % A.LineColor   = [0 0.4470 0.7410];
% % A.violinColor = [0 0.4470 0.7410];
% % A.violinAlpha = 0.4;
% ylabel('$\Delta v$ [mm/s]')
% grid on
% box on
% xticklabels({'Rec 1-imp','fmin 1-imp','Rec 2-imp','fmin 2-imp','Rec J2','fmin J2'})
%% Define colors
col1 = [0.9290 0.6940 0.1250];
col2 = [0.4940 0.1840 0.5560];
col3 = [0.4660 0.6740 0.1880];
col4 = [0 0.4470 0.7410];
colorsRec = [linspace(col3(1),col4(1),6)', linspace(col3(2),col4(2),6)', linspace(col3(3),col4(3),6)'];
colorsNlp = [linspace(col1(1),col2(1),6)', linspace(col1(2),col2(2),6)', linspace(col1(3),col2(3),6)'];
colors = [colorsRec;colorsNlp];


%% Scatter plots
% colororder(colors)
% scatter(1:n,log10(PoCRec),3,'filled')
% hold on
% if alsoNlp
%     scatter(1:n,log10(PoCNlp),3,'filled')
% end
% axis tight
% xlabel('Conjunction ID [-]')
% ylabel('Log of PoC after maneuver [-]')
% grid on
% hold off
% 
% figure()
% if alsoNlp; subplot(1,2,1); end
% colororder(colors(1:6,:))
% a = scatter(1:n,simTimeRec,4,'filled');
% axis tight
% xlabel('Conjunction ID [-]')
% ylabel('Computation time [s]')
% grid on
% ylim([0.15,3.5])
% legend('$n=2$','$n=3$','$n=4$','$n=5$','$n=6$','$n=7$','Interpreter','Latex','color','none','edgeColor','none')
% box on
% 
% if alsoNlp
%     subplot(1,2,2);
%     scatter(1:n,simTimeNlp,4,'filled')
%     axis tight
%     xlabel('Conjunction ID [-]')
%     grid on
%     ylim([0.15,3.5])
%     set(gca,'ColorOrder',colors(7:12,:))
%     legend('$n=2$','$n=3$','$n=4$','$n=5$','$n=6$','$n=7$','Interpreter','Latex','color','none','edgeColor','none')
%     box on
%     yticklabels([])
% end
% 
% figure()
% colororder(colors)
% scatter(1:n,dvRec,4,'filled')
% hold on
% scatter(1:n,dvNlp,4,'filled')
% axis tight
% xlabel('Conjunction ID [-]')
% ylabel('$||\Delta v||$ [mm/s]')
% grid on

edgesInP = -8:1.6:8;
% Build edges
for kk = 1:6
    if alsoNlp
        histinPNlp(kk,:) = histcounts(inPlaneNlp(kk,:),edgesInP);
    end
    histinPRec(kk,:) = histcounts(inPlaneRec(kk,:),edgesInP);
end
if alsoNlp
    histinPNlp  = histinPNlp'/n*100;
else; histinPNlp = []; 
end
histinPRec  = histinPRec'/n*100;

edgesOutP = -0.1:.02:0.1;
% Build edges
for kk = 1:6
    if alsoNlp
        histoutPNlp(kk,:) = histcounts(outPlaneNlp(kk,:),edgesOutP);
    end
    histoutPRec(kk,:) = histcounts(outPlaneRec(kk,:),edgesOutP);
end
if alsoNlp
    histoutPNlp  = histoutPNlp'/n*100;
else; histoutPNlp = []; 
end
histoutPRec  = histoutPRec'/n*100;

edgesNorm = 0:20:200;
% Build edges
for kk = 1:6
    if alsoNlp
        histNormNlp(kk,:) = histcounts(normNlp(kk,:),edgesNorm);
    end
    histNormRec(kk,:) = histcounts(normRec(kk,:),edgesNorm);
end
if alsoNlp
    histNormNlp  = histNormNlp'/n*100;
else; histNormNlp = []; 
end
histNormRec  = histNormRec'/n*100;

% Plot histograms
figure
colororder([colors(1,:);colors(end,:)])
subplot(1,3,1)
bar(edgesNorm(2:end),[histNormRec(:,5),histNormNlp(:,5)],'FaceColor','flat','EdgeColor','flat');
grid on
xticks([20,100,200])
ylim([0,63])
xlabel('$\Delta v$ [mm/s]')
ylabel('\% of cases [-]')
legend('Recursive','fmincon','interpreter','latex','box','off','Orientation','horizontal')
subplot(1,3,2)
bar(edgesInP(2:end),[histinPRec(:,5),histinPRec(:,5)],'FaceColor','flat','EdgeColor','flat');
grid on
yticklabels([])
ylim([0,63])
xlabel('In-plane [deg]')
subplot(1,3,3)
bar(edgesOutP(2:end),[histoutPRec(:,5),histoutPNlp(:,5)],'FaceColor','flat','EdgeColor','flat');
yticks(0:10:60)
yticklabels([])
ylim([0,63])
grid on
xlabel('Out-of-plane [deg]')
