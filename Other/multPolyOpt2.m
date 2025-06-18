initializePath();

%% Initialization variables
load("./data/dataConjunctionsESA.mat");
dvs1  = nan(3,500);
dvs2  = nan(3,500);
poc  = nan(500,1);
md  = nan(500,1);    
compTime  = nan(500,1);    
tcaNewDelta  = nan(500,1);    
smaError = nan(500,1);    
man1 = (6:6:360)/360;
nn = length(man1);
for indCase = 707:2170
for k = 1:nn
    k
    t_man = [man1(k)];
    multiple = 0;
    returnTime = -1;                                                                 % [-] or [days] (1,N) in orbit periods if Earth orbit, days if cislunar
    pp = initOpt(0,0,indCase);
    pp.cislunar = 0;
    pp = defineParams(pp,t_man,returnTime,1e-6,1.5,0);    
    N  = pp.N;                                                                      % [-] (1,1) Number of nodes in the propagation
    n_man = pp.n_man;                                                               % [-] (1,1) Number of nodes where the maneuver can be performed
    if N == 1 && pp.lowThrust; error(['The algorithm needs ' ...
            'at least two nodes to define the low-thrust window']); end
    
    m     = 3  - 2*pp.fixedDir;  pp.m  = m;                                         % [-] (1,1)  Number of optimization variables per node
    u     = zeros(m,n_man);                                                         % [-] (m,N)  Ctrl of the unperturbed trajectory
    
    %% Propagation
    tic;
    [lim,coeff,timeSubtr,xBall,xRetBall] = propDA(pp.DAorder,u,0,pp);
    grad(:,k,indCase) = coeff.C(2:4);

    % %% Optimization
    % if strcmpi(pp.solvingMethod,'lagrange')
    %         [yF,iters] = computeCtrlActiveSet(coeff,u,pp);
    %         % [yF,iters,er,Ys] = computeCtrlRecursive(coeff,u,pp);
    % elseif strcmpi(pp.solvingMethod,'fmincon')
    %         yF = computeCtrlNlp(coeff,u,pp);
    % else
    %     error('The solving method should be either lagrange, or fmincon')
    % end
    % 
    % if pp.fixedDir                                                                  % [-] (3,n_man) Build control matrix node-wise in the case of fixed direction
    %     ctrl = yF.*pp.thrustDirections(:,1:n_man);
    % else
    %     ctrl = yF;                                                                  % [-] (3,n_man) Build control matrix node-wise in the general case
    % end
    % simTime = toc - timeSubtr;     
    % ctrl = pp.ctrlMax*ctrl;
    % 
    % %% Validation
    % metricValPoly = eval_poly(coeff(1).C,coeff(1).E,reshape(yF,1,[]), ...    
    %                             pp.DAorder);
    % [~,~,~,x,xRetMan,x_sec,deltaTca] = propDA(1,ctrl,1,pp);                      % Validate the solution by forward propagating and computing the real PoC
    % if pp.pocType ~= 3
    %     metricValPoly = 10^metricValPoly;
    %     lim           = 10^lim;
    % end

    % lim           = 10^lim;
    % dvs1(:,k,indCase) = ctrl(:,1)*pp.Vsc*1e6;
    % dvs2(:,indCase) = ctrl(:,2)*pp.Vsc*1e6;
    % STMp   = CWStateTransition(pp.primary.n^(3/2),deltaTca/pp.Tsc,0,1);
    % STMs   = CWStateTransition(pp.secondary.n^(3/2),deltaTca/pp.Tsc,0,1);
    % Cpprop = STMp*pp.Cp*STMp';
    % Csprop = STMs*pp.Cs*STMs';
    % Pp     = Cpprop(1:3,1:3);
    % Ps     = Csprop(1:3,1:3);
    % r2ep   = rtn2eci(x(1:3),x(4:6));
    % r2es   = rtn2eci(x_sec(1:3),x_sec(4:6));
    % P      = r2ep*Pp*r2ep' + r2es*Ps*r2es';
    % e2b    = eci2Bplane(x(4:6),x_sec(4:6));
    % e2b    = e2b([1 3],:);
    % PB(:,:,indCase)     = e2b*P*e2b';
    % p      = e2b*(x(1:3)-x_sec(1:3));
    % smd    = dot(p,PB(:,:,indCase)\p);
    % md(k,indCase)     = norm(p)*pp.Lsc;
    % PoC(j) = maximumPc(p,PB(:,:,j),pp.HBR);                                        % [-] (1,1) PoC computed with Chan's formula
    % poc(indCase) = poc_Chan(pp.HBR,PB(:,:,indCase),smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
    % compTime(indCase) = simTime;
    % tcaNewDelta(indCase) = deltaTca;
    % finalCoeMan   = cartesian2kepler(xRetMan,1);     
    % finalCoeBall  = cartesian2kepler(xRetBall,1);
    % finalMeanCoeMan  = osculating2mean(finalCoeMan,1,pp.Lsc);
    % finalMeanCoeBall = osculating2mean(finalCoeBall,1,pp.Lsc);
    % smaError(indCase) = (finalMeanCoeMan.a-finalMeanCoeBall.a)*pp.Lsc*1e3;
end
end
save('gradPoC')
%%
clearvars mus sigmas
% gg(1,:,:) = grad(2,:,:);
% gg(2,:,:) = grad(3,:,:);
% gg(3,:,:) = grad(1,:,:);
gg = grad;
for k = 1:indCase-1
    gN(:,k) = normOfVec(squeeze(gg(:,:,k)));
    gn(:,:,k) = gg(:,:,k);
    gnN(:,k) = gN(:,k)/max(gN(:,k));
end

for i = 1:60
    n = fitdist(gnN(i,:)','normal');
    musN(i) = n.mu;
    sigmasN(i) = n.sigma;
end
x = linspace(6/360,1,60);
% figure()
% plot(x,musN+sigmasN,'k--','HandleVisibility','off')
% hold on
% plot(x,musN-sigmasN,'k--','HandleVisibility','off')
% patch([x flip(x)], [musN-sigmasN flip(musN+sigmasN)],[0.7,0.7,0.7],'DisplayName','1-sigma')
% plot(x,musN,'k','DisplayName','mean value')
% legend('Location','northwest')
% hold off
% xlabel('Orbits before TCA')
% ylabel('Normalized MD gradient')
% axis tight
% grid on
% box on

%%
% col = [0, 0.447, 0.741];
col = 0*ones(1,3);
for k = 1:indCase-1
    gn(:,:,k) = gg(:,:,k)/max(gN(:,k));
end

for j = 1:3
    for i = 1:60
    n(j) = fitdist(squeeze(abs(gn(j,i,:))),'normal');
    mus(j,i) = n(j).mu;
    sigmas(j,i) = n(j).sigma;
    end
end
x = linspace(6/360,1,60);

figure()
lab = ['R','T','N'];

for j = 1:3
    subplot(4,1,j)
    patch([x flip(x)], [max(mus(j,:)-sigmas(j,:),0) flip(mus(j,:)+sigmas(j,:))],col,'FaceAlpha',0.1,'DisplayName','1-$\sigma$','EdgeColor',col)
    hold on
    plot(x,mus(j,:),'k','DisplayName','Mean value','Color',col)
    hold off
    ylabel(lab(j))
    axis tight
    grid on
    box on
    xticklabels('')
    ylim([0,1.13])
end
subplot(4,1,4)
patch([x flip(x)], [max(musN-sigmasN,0) flip(musN+sigmasN)],col,'FaceAlpha',0.1,'DisplayName','1-$\sigma$','EdgeColor',col)
hold on
plot(x,musN,'k','DisplayName','mean value','Color',col)
legend('Location','northwest')
hold off
xlabel('Orbits before TCA')
ylabel('Norm')
axis tight
grid on
box on
xlabel('Orbits before TCA')
legend('Location','northwest','Interpreter','latex')

%%
% close all
% sz = 10;
% figure
% subplot(2,1,1)
% scatter(man1,dvs1'/1e3,sz,'filled')
% hold on
% scatter(man1,normOfVec(dvs1)/1e3,sz,[0 0 0],'filled')
% hold off
% ylabel('$\Delta v_1$ [m/s]')
% legend('R','T','N','$|\cdot|$','interpreter','latex')
% set(gca, 'XDir', 'reverse')
% box on
% grid on
% axis tight
% xticks(0:0.5:5.5);
% xticklabels('')
% 
% subplot(2,1,2)
% scatter(man1,dvs2'/1e3,sz,'filled')
% hold on
% scatter(man1,normOfVec(dvs2)/1e3,sz,[0 0 0],'filled')
% hold off
% xlabel('Number of orbits from TCA [-]')
% ylabel('$\Delta v_2$ [m/s]')
% set(gca, 'XDir', 'reverse')
% box on
% grid on
% axis tight
% xticks(0:0.5:5.5);
% 
% figure
% subplot(2,1,1)
% scatter(man1,md,sz,'filled')
% ylabel('MD [km]')
% box on
% grid on
% set(gca, 'XDir', 'reverse')
% axis tight
% xticklabels('')
% xticks(0:0.5:5.5);
% 
% subplot(2,1,2)
% scatter(man1,log10(abs(smaError)),sz,'filled')
% xlabel('Number of orbits from TCA [-]')
% ylabel('SMA error [m]')
% box on
% grid on
% set(gca, 'XDir', 'reverse')
% yticks(-1:3)
% yticklabels({'0.1','1','10','100','1000'})
% axis tight
% xticks(0:0.5:5.5);
