%% load variables
close all
n = 2170;
kk = 5;
load('SimOutput\singleImpulseEsaCases\convex.mat')
PoCConv     = PoC;
simTimeConv = simTime;
dvConv = dv*1000;
load(['SimOutput\gridRefinedEsaCases\nlp',num2str(kk),'.mat'])
PoCNlp(1,:)     = PoC;
simTimeNlp(1,:) = simTime;
dvNlp(1,:)      = normOfVec(dvs);
% load(['SimOutput\gridRefinedEsaCases\nlp',num2str(kk),'_2Imp.mat'])
% PoCNlp(2,:)     = PoC;
% simTimeNlp(2,:) = simTime;
% dvNlp(2,:)      = normOfVec(dvs);
load(['SimOutput\gridRefinedEsaCases\recursive',num2str(kk),'.mat'])
PoCRec(1,:)     = PoC;
simTimeRec(1,:) = simTime;
dvRec(1,:)      = normOfVec(dvs);
% load(['SimOutput\gridRefinedEsaCases\recursive',num2str(kk),'_2Imp.mat'])
% PoCRec(2,:)     = PoC;
% simTimeRec(2,:) = simTime;
% dvRec(2,:)      = normOfVec(dvs);
clearvars -except PoCConv simTimeConv PoCNlp simTimeNlp PoCRec simTimeRec dvNlp dvRec dvConv nodeThrust n

%% Define colors
col1      = [0.9290 0.6940 0.1250];
col2      = [0.4940 0.1840 0.5560];
col3      = [0.4660 0.6740 0.1880];
col4      = [0 0.4470 0.7410];
colorsRec = [linspace(col3(1),col4(1),6)', linspace(col3(2),col4(2),6)', linspace(col3(3),col4(3),6)'];
colorsNlp = [linspace(col1(1),col2(1),6)', linspace(col1(2),col2(2),6)', linspace(col1(3),col2(3),6)'];
colors    = [colorsRec([1,end],:);colorsNlp([1,end],:)];

%% Scatter plots
figure()
colororder(colors)
PoCRec(PoCRec<1e-8) = nan;
PoCNlp(PoCNlp<1e-8) = nan;
PoCConv(PoCConv<1e-8) = nan;
scatter(1:n,log10(PoCRec),3,'filled')
hold on
scatter(1:n,log10(PoCNlp),3,'filled')
scatter(1:n,log10(PoCConv),3,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('PoC after maneuver [-]')
grid on
hold off
legend('Recursive','fmincon')

figure()
colororder(colors)
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
colororder(colors)
scatter(1:n,dvRec,4,'filled')
hold on
scatter(1:n,dvNlp,4,'filled')
axis tight
xlabel('Conjunction ID [-]')
ylabel('\Delta v [mm/s]')
grid on
legend('Recursive','fmincon')

figure()
colororder(colors)
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

%% Spiral for thrusting times
figure()
t = linspace(0,2.5,5000)*2*pi;
tt = linspace(0.5,2.5,11)*2*pi;
for i = 1:length(nodeThrust)
    piThrust(i) = tt(nodeThrust(i));
end
for i = 1:length(tt)
    a(i)     = sum(piThrust==tt(i));
    [~,b(i)] = (min(abs(tt(i)-t)));
end
r = 1; 
x = r*cos(t);
y = r*sin(t);
z = t/2/pi;
plot3(x,y,z,'k');
hold on
for i = 1:length(tt)
    plot3(x(b(i)),y(b(i)),z(b(i)),'ko')
end
for i = 1:length(tt)
    if a(i) > 0
        plot3(x(b(i)),y(b(i)),z(b(i)),'.','color',[0 0 0],'MarkerSize',a(i)/max(a)*100)
    end
end
plot3(x(1),y(1),z(1),'ro')
zlabel('Number of orbits [-]')
xticks([])
yticks([])
grid on