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
load('SimOutput\IAC\rec5.mat')
PoC(PoC>prctile(PoC,99)) = nan;
PoC(PoC<prctile(PoC,1)) = nan;
PoCRec     = PoC;
simTimeRec = compTime;
errPRec = errP;
errVRec = errV;
dv1     = squeeze(dvs(:,1,:));
dv2     = squeeze(dvs(:,2,:));
dv3     = squeeze(dvs(:,3,:));
dvn1 = normOfVec(dv1);
dvn2 = normOfVec(dv2);
dvn3 = normOfVec(dv3);
% dvn(dvn>prctile(dvn,95)) = nan;
% inPlaneRec(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
% outPlaneRec(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';

%% Violin plots
figure
A = iosr.statistics.boxPlot(PoCRec');
A.mediancolor = 'b';
A.outliersize = nan;
A.showMean = true;
A.percentile = [50 50]; 
A.lineStyle = 'none';
A.meancolor = 'b';
A.showViolin = true;
A.showScatter = false;
A.LineColor   = [0 0.4470 0.7410];
A.violinColor = [0 0.4470 0.7410];
A.violinAlpha = 0.4;
xlabel('Expansion order [-]')
ylabel('PoC [-]')
grid on
box on

%% Scatter plots
figure()
scatter(1:n,log10(PoCRec),4,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('Log of PoC after maneuver [-]')
grid on
hold off

figure()
a = scatter(1:n,simTimeRec,4,'filled');
axis tight
xlabel('Conjunction ID [-]')
ylabel('Computation time [s]')
grid on
ylim([2,5])
box on

edgesNorm = 0:20:1000;
% Build edges
histNormRec(1,:) = histcounts(dvn1,edgesNorm);
histNormRec(2,:) = histcounts(dvn2,edgesNorm);
histNormRec(3,:) = histcounts(dvn3,edgesNorm);

histNormRec  = histNormRec'/n*100;

% Plot histograms
figure
bar(edgesNorm(2:end),histNormRec','FaceColor','flat','EdgeColor','flat');
grid on
xlabel('$\Delta v$ [mm/s]')
ylabel('\% of cases [-]')
legend('$\Delta v_1$','$\Delta v_2$','$\Delta v_3$','interpreter','latex','box','off','Orientation','horizontal')

