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
close all
n = 2170;
PoCNlp      = nan(6,n);
PoCRec      = nan(6,n);
simTimeNlp  = nan(6,n);
simTimeRec  = nan(6,n);
dvNlp       = nan(6,n);
dvRec       = nan(6,n);
load('SimOutput\singleImpulseEsaCases\convex.mat')
PoCConv     = PoC;
simTimeConv = simTime;
for kk = 2:6
    load(['SimOutput\doubleImpulseEsaCases\nlp',num2str(kk),'.mat'])
    PoC(PoC>prctile(PoC,99)) = nan;
    PoC(PoC<prctile(PoC,1)) = nan;
    PoCNlp(kk-1,:)     = PoC;
    simTimeNlp(kk-1,:) = simTime;
    dvNlp(kk-1,:)      = normOfVec(dvs);
    load(['SimOutput\doubleImpulseEsaCases\rec',num2str(kk),'.mat'])
    PoC(PoC>prctile(PoC,99)) = nan;
    PoC(PoC<prctile(PoC,1)) = nan;
    PoCRec(kk-1,:)     = PoC;
    simTimeRec(kk-1,:) = simTime;
    dvRec(kk-1,:)      = normOfVec(dvs);
end
clearvars -except PoCConv simTimeConv PoCNlp simTimeNlp PoCRec simTimeRec dvNlp dvRec n
alsoNlp = 1;
%% Violin plots
% figure
% A = iosr.statistics.boxPlot(PoCRec');
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
%     B = iosr.statistics.boxPlot(PoCNlp');
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
% figure
% subplot(1,2,1)
% C = iosr.statistics.boxPlot(2,PoCRec(1,:)');
% C.mediancolor = 'b';
% C.outliersize = nan;
% C.showMean = true;
% C.percentile = [50,50]; 
% C.lineStyle = 'none';
% C.meancolor = 'b';
% C.showViolin = true;
% C.showScatter = false;
% C.LineColor   = [0 0.4470 0.7410];
% C.violinColor = [0 0.4470 0.7410];
% C.violinAlpha = 0.4;
% if alsoNlp 
%     hold on
%     D = iosr.statistics.boxPlot(2,PoCNlp(1,:)');
%     D.mediancolor = 'r';
%     D.outliersize = nan;
%     D.showMean = true;
%     D.percentile = [50,50]; 
%     D.lineStyle = 'none';
%     D.meancolor = 'r';
%     D.showViolin = true;
%     D.showScatter = false;
%     D.LineColor = [0.4940 0.1840 0.5560];
%     D.violinColor = [0.4940 0.1840 0.5560];
%     D.violinAlpha = 0.3;
% end
% xlabel('Expansion order [-]')
% ylabel('PoC [-]')
% grid on
% box on
% subplot(1,2,2)
% C = iosr.statistics.boxPlot(3,PoCRec(2,:)');
% C.mediancolor = 'b';
% C.outliersize = nan;
% C.showMean = true;
% C.percentile = [50,50]; 
% C.lineStyle = 'none';
% C.meancolor = 'b';
% C.showViolin = true;
% C.showScatter = false;
% C.LineColor   = [0 0.4470 0.7410];
% C.violinColor = [0 0.4470 0.7410];
% C.violinAlpha = 0.4;
% if alsoNlp
%     hold on
%     D = iosr.statistics.boxPlot(3,PoCNlp(2,:)');
%     D.mediancolor = 'r';
%     D.outliersize = nan;
%     D.showMean = true;
%     D.percentile = [50,50]; 
%     D.lineStyle = 'none';
%     D.meancolor = 'r';
%     D.showViolin = true;
%     D.showScatter = false;
%     D.LineColor = [0.4940 0.1840 0.5560];
%     D.violinColor = [0.4940 0.1840 0.5560];
%     D.violinAlpha = 0.3;
% end
% grid on
% box on
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

figure()
if alsoNlp; subplot(1,2,1); end
colororder(colors(1:6,:))
a = scatter(1:n,simTimeRec,4,'filled');
axis tight
xlabel('Conjunction ID [-]')
ylabel('Computation time [s]')
grid on
ylim([0.15,2])
legend('$n=2$','$n=3$','$n=4$','$n=5$','$n=6$','$n=7$','Interpreter','Latex')
box on

if alsoNlp
    subplot(1,2,2);
    scatter(1:n,simTimeNlp,4,'filled')
    axis tight
    xlabel('Conjunction ID [-]')
    grid on
    ylim([0.15,2])
    set(gca,'ColorOrder',colors(7:12,:))
    legend('$n=2$','$n=3$','$n=4$','$n=5$','$n=6$','$n=7$','Interpreter','Latex')
    box on
    yticklabels([])
end

figure()
colororder(colors)
scatter(1:n,dvRec,4,'filled')
hold on
scatter(1:n,dvNlp,4,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('$||\Delta v||$ [mm/s]')
grid on

%% Histograms
% Define edges
edgesTime = 0.1:0.05:0.45;
edgesPoC = -6.0002:0.00002:-5.9998;
edgesDv = 5:10:1000;
% Build edges
for kk = 1:6
    if alsoNlp
        histTimeNlp(kk,:) = histcounts(simTimeNlp(kk,:),edgesTime);
        histPocNlp(kk,:) = histcounts(log10(PoCNlp(kk,:)),edgesPoC);
        histDvNlp(kk,:) = histcounts(dvNlp(kk,:),edgesDv);
    end
    histTimeRec(kk,:) = histcounts(simTimeRec(kk,:),edgesTime);
    histPocRec(kk,:) = histcounts(log10(PoCRec(kk,:)),edgesPoC);
    histDvRec(kk,:) = histcounts(dvRec(kk,:),edgesDv);
end
histPocConv  = histcounts(log10(PoCConv),edgesPoC);
histPocConv  = histPocConv'/n*100;
if alsoNlp
    histTimeNlp = histTimeNlp'/n*100;
    histPocNlp  = histPocNlp'/n*100;
    histDvNlp  = histDvNlp'/n*100;
else; histTimeNlp = []; histPocNlp = []; histDvNlp = []; 
end
histTimeRec = histTimeRec'/n*100;
histPocRec  = histPocRec'/n*100;
histDvRec  = histDvRec'/n*100;

% Plot histograms
% figure
% colororder(colors)
% bar(edgesTime(2:end),[histTimeRec,histTimeNlp],'FaceColor','flat');
% legend
% grid on
% xlabel('Computation time [s]')
% ylabel('\% of cases [-]')
% legend('$n=2$','$n=3$','$n=4$','$n=5$', ...
%        '$n=6$','$n=7$','$n=2$','$n=3$', ...
%        '$n=4$','$n=5$', '$n=6$','$n=7$','interpreter','latex')

% figure
% colororder(colors)
% bar(edgesPoC(2:end),[histPocRec,histPocNlp,histPocConv],'grouped','FaceColor','flat');
% grid on
% xlabel('log$_{10}$(PoC) [-]')
% ylabel('\% of cases [-]')
% legend('$n=2$','$n=3$','$n=4$','$n=5$', ...
%        '$n=6$','$n=7$','$n=2$','$n=3$', ...
%        '$n=4$','$n=5$', '$n=6$','$n=7$','interpreter','latex')

figure
colororder(colors)
bar(edgesDv(2:end),[histDvRec,histDvNlp],'FaceColor','flat');
legend
grid on
xlabel('$\Delta v$ [mm/s]')
ylabel('\% of cases [-]')
legend('$n=2$','$n=3$','$n=4$','$n=5$', ...
       '$n=6$','$n=7$','$n=2$','$n=3$', ...
       '$n=4$','$n=5$', '$n=6$','$n=7$','interpreter','latex','box','off')

% Plot normal distribution