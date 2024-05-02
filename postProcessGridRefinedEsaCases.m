%% load variables
% close all
n = 2170;
kk = 5;
load('SimOutput\singleImpulseEsaCases\convex.mat')
PoCConv     = PoC;
simTimeConv = simTime;
load(['SimOutput\gridRefinedEsaCases\nlp',num2str(kk),'_2Imp.mat'])
PoCNlp     = PoC;
simTimeNlp = simTime;
% dvNlp      = normOfVec(dvs)*sqrt(pp.mu/6378)*1e6;
dvNlp      = normOfVec(dvs);
load(['SimOutput\gridRefinedEsaCases\recursive',num2str(kk),'_2Imp.mat'])
PoCRec     = PoC;
simTimeRec = simTime;
% dvRec      = normOfVec(dvs)*sqrt(pp.mu/6378)*1e6;
dvRec      = normOfVec(dvs);
clearvars -except PoCConv simTimeConv PoCNlp simTimeNlp PoCRec simTimeRec dvNlp dvRec nodeThrust n


%% Scatter plots
figure()
PoCRec(PoCRec<1e-8) = nan;
PoCNlp(PoCNlp<1e-8) = nan;
scatter(1:n,log10(PoCRec),3,'filled')
hold on
scatter(1:n,log10(PoCNlp),3,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('PoC after maneuver [-]')
grid on
hold off
legend('Recursive','fmincon')

figure()
scatter(1:n,simTimeRec,4,'filled')
hold on
scatter(1:n,simTimeNlp,4,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('Computation time [-]')
grid on
ylim([0.3,1])
legend('Recursive','fmincon')

figure()
scatter(1:n,nodeThrust,3,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('PoC after maneuver [-]')
grid on
hold off

%% Histograms
% Define edges
edgesTime = 0.1:0.05:0.45;
edgesPoC = -6.04:0.01:-5.96;
edgesDv = (0:3e-6:5e-5)*1e6*sqrt(398600/6378);
% Build edges
histTimeNlp = histcounts(simTimeNlp,edgesTime);
histTimeRec = histcounts(simTimeRec,edgesTime);
histPocNlp  = histcounts(log10(PoCNlp),edgesPoC);
histPocRec  = histcounts(log10(PoCRec),edgesPoC);
histDvNlp   = histcounts(dvNlp,edgesDv);
histDvRec   = histcounts(dvRec,edgesDv);
histTimeNlp = histTimeNlp'/n*100;
histTimeRec = histTimeRec'/n*100;
histPocNlp  = histPocNlp'/n*100;
histPocRec  = histPocRec'/n*100;
histDvNlp  = histDvNlp'/n*100;
histDvRec  = histDvRec'/n*100;

% Plot histograms
figure
bar(edgesTime(2:end),[histTimeRec,histTimeNlp],'FaceColor','flat');
legend
grid on
xlabel('Computation time [s]')
ylabel('\% of cases [-]')

figure
bar(edgesPoC(2:end),[histPocRec,histPocNlp],'grouped','FaceColor','flat');
grid on
xlabel('log$_{10}$(PoC) [-]')
ylabel('\% of cases [-]')

figure
bar(edgesDv(2:end),[histDvRec,histDvNlp],'FaceColor','flat');
legend
grid on
xlabel('$\Delta v$ [mm/s]')
ylabel('\% of cases [-]')

