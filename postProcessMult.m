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
y = 2:7;
yy = length(y);
% PoCNlp      = nan(yy,n);
% PoCRec      = nan(yy,n);
% simTimeNlp  = nan(yy,n);
% simTimeRec  = nan(yy,n);
% normNlp     = nan(yy,n);
% normRec     = nan(yy,n);
% inPlaneNlp  = nan(yy,n);
% inPlaneRec  = nan(yy,n);
% outPlaneNlp = nan(yy,n);
% outPlaneRec = nan(yy,n);
% for kk = y
%     % PoC(PoC>prctile(PoC,95)) = nan;
%     % PoC(PoC<prctile(PoC,5)) = nan;
%     % % tcaNewDelta(tcaNewDelta>prctile(tcaNewDelta,99.95)) = nan;
%     % % tcaNewDelta(tcaNewDelta<prctile(tcaNewDelta,0.05)) = nan;
%     % % compTime(compTime>prctile(compTime,99)) = nan;
%     % PoCNlp(kk-1,:)     = PoC;
%     % tcaRec(kk-1,:)     = nan(2170,1);
%     % compTime(compTime>prctile(compTime,99)) = nan;
%     % simTimeNlp(kk-1,:) = compTime;
%     % dvs = squeeze(dvs);
%     % inPlaneNlp(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
%     % outPlaneNlp(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
%     % dvn1 = normOfVec(squeeze(dvs));
%     % dvn1(dvn1>prctile(dvn1,99.9)) = nan;
%     % dv1 = dvs;
%     % normNlp(kk-1,:)     = dvn1;
%     load('SimOutput\IACReturn.mat')
%     PoC(PoC>prctile(PoC,95)) = nan;
%     PoC(PoC<prctile(PoC,5)) = nan;
%     tcaNewDelta(tcaNewDelta>prctile(tcaNewDelta,99.95)) = nan;
%     tcaNewDelta(tcaNewDelta<prctile(tcaNewDelta,0.05)) = nan;
%     tcaNlp(kk-1,:)     = tcaNewDelta;
%     compTime(compTime>prctile(compTime,99)) = nan;
%     PoCRec(kk-1,:)     = PoC;
%     simTimeRec(kk-1,:) = compTime;
%     % dvs = squeeze(dvs);
%     % inPlaneRec(kk-1,:)  = rad2deg(atan(dvs(1,:)'./dvs(2,:)'))';
%     % outPlaneRec(kk-1,:) = rad2deg(atan(dvs(3,:)'./normOfVec(dvs(1:2,:))'))';
%     % dvn = normOfVec(squeeze(dvs));
%     % dvn(dvn>prctile(dvn,99.9)) = nan;
%     % normRec(kk-1,:)     = dvn;  
%     % cr(:,:,kk-1)     = convRad*pp.scaling(4)*pp.ctrlMax*1e6;  
% end
load('SimOutput\IACReturn');
for j = 1:2170
    dv1(j) = sum(normOfVec(dvs(:,:,j)));
end
valid     = boolean(sum(iterationsN,1)<pp.maxIter);
dv1       = dv1(valid);
PoC1      = PoC(valid);
compTime1 = compTime(valid);
tcaRec1   = tcaNewDelta(valid);
meanAErr1 = meanAErr(valid);
meanEErr1 = meanEErr(valid);
rRetErr1  = rRetErr(valid);
vRetErr1  = vRetErr(valid);
dvTot     = vRetErr(valid);

load('SimOutput\IACMean');
for j = 1:2170
    dv2(j) = sum(normOfVec(dvs(:,:,j)));
end
valid     = boolean(meanAErr<5);
dv2       = dv2(valid);
PoC2      = PoC(valid);
compTime2 = compTime(valid);
tcaRec2   = tcaNewDelta(valid);
meanAErr2 = meanAErr(valid);
meanEErr2 = meanEErr(valid);
rRetErr2  = rRetErr(valid);
vRetErr2  = vRetErr(valid);
alsoNlp   = 0;

%% Violin plots
% figure
% A = violin([PoC1', PoC2']*1e6-1);
% hold off
% xlabel('Expansion order [-]')
% ylabel('PoC relative error [-]')

% figure
% subplot(3,1,1)
% A = violin([compTime1', compTime2']);
% hold off
% xlabel('Expansion order [-]')
% ylabel('Computation Time [s]')


%% Define colors
% col1 = [0.9290 0.6940 0.1250];
% col2 = [0.4940 0.1840 0.5560];
% col3 = [0.4660 0.6740 0.1880];
% col4 = [0 0.4470 0.7410];
% colorsRec = [linspace(col3(1),col4(1),yy)', linspace(col3(2),col4(2),yy)', linspace(col3(3),col4(3),yy)'];
% colorsNlp = [linspace(col1(1),col2(1),yy)', linspace(col1(2),col2(2),yy)', linspace(col1(3),col2(3),yy)'];
% colors = [colorsRec;colorsNlp];
% colororder(colors)

% crR = squeeze(cr(1,:,:));
% crT = squeeze(cr(2,:,:));
% crN = squeeze(cr(3,:,:));

% figure()
% [crRsorted,ord] = sort(crR);
% semilogy(abs(dvs(1,ord(:,4))),'.','LineWidth',1)
% hold on
% semilogy(crRsorted,'LineWidth',2)
% axis tight
% box on
% grid on
% hold off
% ylabel('$\Delta v_R$ [mm/s]')
% xlabel('Conjunction ID [-]')
%     % 
    % figure()
    % [crTsorted,ord] = sort(crT);
    % ddv = abs(dvs(2,ord(:,4)));
    % upB = crTsorted(:,5)';
    % ddvup = ddv; ddvup(ddvup<upB*2) = nan;
    % semilogy(abs(ddv),'b.','LineWidth',1)
    % hold on
    % semilogy(abs(ddvup),'r.','LineWidth',1)
    % semilogy(crTsorted,'LineWidth',2)
    % legend('','','n=2','n=3','n=4','n=5','n=6','n=7','Interpreter','latex','Orientation','horizontal')
    % axis tight
    % box on
    % grid on
    % hold off
    % ylabel('$\Delta v_T$ [mm/s]')
    % xlabel('Conjunction ID [-]')
    
    % iters = sum(iterationsN(1:2,ord(:,4)),1);
    % iterred = iters; iterred(isnan(ddvup)) = nan;
    % figure()
    % plot(iters,'b.') 
    % hold on 
    % plot(iterred,'r.')
    % hold off
    % axis tight
    % box on
    % grid on
    % hold off
    % ylabel('Iterations [-]')
    % xlabel('Conjunction ID [-]')

% figure()
% [crNsorted,ord] = sort(crN);
% semilogy(abs(dvs(3,ord(:,4))),'.','LineWidth',1)
% hold on
% semilogy(crNsorted,'LineWidth',2)
% axis tight
% box on
% grid on
% hold off
% ylabel('$\Delta v_N$ [mm/s]')
% xlabel('Conjunction ID [-]')
% 
% %% Delta V histograms
% figure()
% edgesInP = -8:1.6:8;
% edgesOutP = -0.1:.02:0.1;
% edgesNorm = 0:20:200;
% colororder([colors(1,:);colors(end,:)])
% subplot(1,3,1)
% histograms([normRec(4,:);normNlp(4,:)],edgesNorm)
% xlabel('$\Delta v$ [mm/s]')
% xticks([20,100,200])
% legend('Recursive','fmincon','interpreter','latex','box','off','Orientation','horizontal')
% subplot(1,3,2)
% histograms([inPlaneRec(4,:);inPlaneNlp(4,:)],edgesInP)
% xlabel('In-plane [deg]')
% yticklabels([])
% ylim([0,63])
% subplot(1,3,3)
% histograms([outPlaneRec(4,:);outPlaneNlp(4,:)],edgesOutP)
% xlabel('Out-of-plane [deg]')
% yticks(0:10:60)
% yticklabels([])
% ylim([0,63])

%% TCA histogram
placeFigure
edges = -0.5:0.05:0.5;
histograms(tcaNewDelta,edges)
xlabel('$\Delta t_{CA}$ [s]')


placeFigure
edges = 0.1:0.005:0.3;
histograms(compTime2,edges)
xlabel('Computation time [s]')
legend('Return','Mean')

% 
%% Return histograms
placeFigure
subplot(2,1,1)
edges = 0:0.01:.5;
histograms(rRetErr1,edges)
xlabel('$e_{r}$ [m]')
edges = 0:0.01:.2;
subplot(2,1,2)
histograms(vRetErr1,edges)
xlabel('$e_{v}$ [mm/s]')

placeFigure
subplot(2,1,1)
edges = 0:0.1:5;
histograms(meanAErr2,edges)
xlabel('$e_{a}$ [m]')
subplot(2,1,2)
histograms(meanEErr2,edges)
xlabel('$e_{e}$ [m]')

%% Delta V histograms
% dv1(abs(dv1)>1000) = nan;
% dvDiff = squeeze(dvs) - dv1;
% dv1(:,any(isnan(dv1))) = [];
% dvDiff(abs(dvDiff)>200) = nan;
% dvDiff(:,any(isnan(dvDiff))) = [];
% figure()
% edges =  -.2:.02:.2;
% colororder([colors(1,:);colors(end,:)])
% subplot(1,3,1)
% histograms([dvDiff(1,:)],edges)
% xlabel('R [mm/s]')
% ylim([0,90])
% % legend('R','T','N','interpreter','latex','box','off','Orientation','horizontal')
% subplot(1,3,2)
% edges =  -1:.1:1;
% histograms([dvDiff(2,:)],edges)
% xlabel('T [mm/s]')
% yticklabels([])
% ylim([0,90])
% subplot(1,3,3)
% edges =  -.03:.003:.03;
% histograms([dvDiff(3,:)],edges)
% xlabel('N [mm/s]')
% yticks(0:10:60)
% yticklabels([])
% ylim([0,90])

%%
% dvLog = log10(abs(dvDiff));
% figure
% % A = iosr.statistics.boxPlot(dvLog');
% A = violin(dvLog');
% % A.outliersize = 15;
% % A = violin(log10(abs(dvs))');
% hold on
% % B = violin(log10(abs(dv1))');    
% % B.LineColor = [0.4940 0.1840 0.5560];
% % B.violinColor = [0.4940 0.1840 0.5560];
% % B.mediancolor = 'r';
% % B.meancolor   = 'r';
% yticks([-10:5:0,1])
% % A.percentile = [10 90]; 
% hold off
% xticklabels({'R','T','N'})
% ylabel('$\Delta V$ [log(mm/s)]')

%%