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
for kk = 2:7
    load(['SimOutput\singleImpulseEsaCases\nlp',num2str(kk),'.mat'])
    PoCNlp(kk-1,:)     = PoC;
    simTimeNlp(kk-1,:) = simTime;
    dvNlp(kk-1,:)      = normOfVec(dvs)*sqrt(pp.mu/6378)*1e6;
    load(['SimOutput\singleImpulseEsaCasesFixedDir\rec',num2str(kk),'.mat'])
    PoCRec(kk-1,:)     = PoC;
    simTimeRec(kk-1,:) = simTime;
    dvRec(kk-1,:)      = normOfVec(dvs)*sqrt(pp.mu/6378)*1e6;
end
clearvars -except PoCConv simTimeConv PoCNlp simTimeNlp PoCRec simTimeRec dvNlp dvRec n
alsoNlp = 0;

%% Define colors
col1 = [0.9290 0.6940 0.1250];
col2 = [0.4940 0.1840 0.5560];
col3 = [0.4660 0.6740 0.1880];
col4 = [0 0.4470 0.7410];
colorsRec = [linspace(col3(1),col4(1),6)', linspace(col3(2),col4(2),6)', linspace(col3(3),col4(3),6)'];
colorsNlp = [linspace(col1(1),col2(1),6)', linspace(col1(2),col2(2),6)', linspace(col1(3),col2(3),6)'];
colors = [colorsRec;colorsNlp];


%% Scatter plots
figure()
PoCRec(PoCRec<1e-8) = nan;
PoCNlp(PoCNlp<1e-8) = nan;
colororder(colors)
scatter(1:n,log10(PoCRec),3,'filled')
hold on
if alsoNlp
    scatter(1:n,log10(PoCNlp),3,'filled')
end
axis tight
xlabel('Conjunction ID [-]')
ylabel('Log of PoC after maneuver [-]')
grid on
hold off

figure()
colororder(colors(1:6,:))
scatter(1:n,simTimeRec,4,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('Computation time [s]')
grid on
ylim([0.15,0.4])

if alsoNlp
    figure()
    colororder(colors(7:12,:))
    scatter(1:n,simTimeNlp,4,'filled')
    axis tight
    xlabel('Conjunction ID [-]')
    ylabel('Computation time [s]')
    grid on
    ylim([0.15,0.4])
end

%% Histograms
% Define edges
edgesTime = 0.1:0.05:0.45;
edgesPoC = -6.0002:0.00002:-5.9998;
edgesDv = (0:3e-6:5e-5)*1e6*sqrt(398600/6378);
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
if alsoNlp
    histTimeNlp = histTimeNlp'/n*100;
    histPocNlp  = histPocNlp'/n*100;
    histDvNlp  = histDvNlp'/n*100;
else; histTimeNlp = []; histPocNlp = []; histDvNlp = []; end
histTimeRec = histTimeRec'/n*100;
histPocRec  = histPocRec'/n*100;
histDvRec  = histDvRec'/n*100;

% Plot histograms
figure
colororder(colors)
bar(edgesTime(2:end),[histTimeRec,histTimeNlp],'FaceColor','flat');
legend
grid on
xlabel('Computation time [s]')
ylabel('\% of cases [-]')

figure
colororder(colors)
bar(edgesPoC(2:end),[histPocRec,histPocNlp],'grouped','FaceColor','flat');
grid on
xlabel('log$_{10}$(PoC) [-]')
ylabel('\% of cases [-]')

figure
colororder(colors)
bar(edgesDv(2:end),[histDvRec,histDvNlp],'FaceColor','flat');
legend
grid on
xlabel('$\Delta v$ [mm/s]')
ylabel('\% of cases [-]')

