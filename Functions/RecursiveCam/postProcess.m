function [] = postProcess(smdLim,xBall,xMan,lim,ctrl,simTime,pp)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
metricFlag = pp.metricFlag;
mdLim      = pp.mdLim;
Lsc        = pp.Lsc;
Vsc        = pp.Vsc;
x_sTCA     = pp.x_sTCA;
P          = pp.P;
%nodes for low-thrust
if pp.lowThrust
    t_lt            = pp.t; 
    t_lt(pp.isConj) = [];
    dt_lt           = diff(t_lt); 
    dt_lt           = abs(dt_lt(pp.canFire));
    s               = strfind(pp.canFire',[1 0]) + 1; % Nodes in which the thruster is shut-off
    if length(s) > 1
        for j = 1:length(s)-1
            tt(:,j) = linspace(pp.t(s(j))-1e-6,pp.t(s(j)+1)+1e-6,2); 
        end
        ttt = [reshape(tt,[],1);pp.t(end-1)-1e-6];
    else
        ttt = pp.t(s)-1e-6;
    end
end
% compute PoC after maneuver
PoC = nan(pp.n_conj,1);
for k = 1:pp.n_conj
    xb     = xBall(:,k);
    x      = xMan(:,k);
    x_s    = x_sTCA(:,k);
    e2b    = eci2Bplane(xb(4:6),x_s(4:6));
    e2b    = e2b([1 3],:);
    PB     = e2b*P(:,:,k)*e2b';
    p      = e2b*(x(1:3)-x_s(1:3));
    smd    = dot(p,PB\p);
    PoC(k) = poc_Chan(pp.HBR(k),PB,smd,3);                                        % [-] (1,1) PoC computed with Chan's formula
end
poc_tot = PoCTot(PoC);
disp(['Solver: ', pp.solvingMethod])
disp(['Computation time ',num2str(simTime), ' s'])
n   = size(ctrl,2);
disp(['Number of conjunctions: ', num2str(pp.n_conj)])
if ~pp.lowThrust
    ctrl  = ctrl*Vsc*1e6;
    dv = ctrl;
    disp(['Number of Delta-vs = ',num2str(n)])
    for i = 1:size(dv,2)
        disp(['Delta-v',num2str(i),' = ',num2str(norm(dv(:,i))), ' mm/s'])
    end
    disp(['Total Delta-v = ',num2str(sum(normOfVec(dv))), ' mm/s'])
else
    disp(['Total Delta-v = ',num2str(normOfVec(ctrl)*dt_lt*Vsc*1e6), ' mm/s'])
end

if metricFlag == 3
    disp(['Limit: ',num2str(sqrt(lim)*Lsc), ' km'])
else
    disp(['Limit: ',num2str(lim)])
end
% disp(['Targeting error position ',num2str(norm(distVal(1:3))*pp.Lsc),' km']);
% disp(['Targeting error velocity ',num2str(norm(distVal(4:6))*pp.Vsc*1e3),' m/s']);
if metricFlag == 0 || metricFlag == 1
    disp(['PoC after validation ',num2str(poc_tot)]);
    stringMetric = 'PoC';
elseif metricFlag == 2
    disp(['SMD after validation ',num2str(metricVal)]);                    % [-] (1,1) SMD limit computed with Alfriend and Akella's formula applied to PC
    stringMetric = 'SMD';
else
    disp(['Miss distance after validation ',num2str(sqrt(metricVal)*Lsc)]);      % [-] (1,1) SMD limit computed with Alfriend and Akella's formula applied to PC
    stringMetric = 'miss distance';
end
disp(['Metric used: ', stringMetric])

ctrlNorm = normOfVec(ctrl);

figure()
if ~pp.lowThrust
    ctrl(ctrl==0) = nan;
    if ~pp.cislunar
        t  = pp.t(pp.canFire)/pp.T;
    else
        t  = pp.t(pp.canFire)*pp.Tsc/86400;
    end
    stem(t,ctrl(1,:),'LineWidth',2)
    hold on
    stem(t,ctrl(2,:),'LineWidth',2)
    stem(t,ctrl(3,:),'LineWidth',2)
    stem(t(ctrlNorm~=0),ctrlNorm(ctrlNorm~=0),'color','k','LineWidth',2)
    plot(t(ctrlNorm==0),ctrlNorm(ctrlNorm==0),'color','k')
    ylabel('$\Delta v$ [mm/s]')
else
    ctrl = ctrl*pp.Asc*1e6;
    for i = 2:pp.N
        if ~pp.canFire(i) && pp.canFire(i-1)
            ctrl = [ctrl(:,1:i-1) zeros(3,1) ctrl(:,i:end)];
        end
    end
    [t,o] = sort([pp.t;ttt],'descend');
    if ~pp.cislunar
        t     = t/pp.T;
    else
        t  = t*pp.Tsc/86400;
    end
    ctrl  = [ctrl, zeros(3,length(ttt)+1)];
    ctrl  = ctrl(:,o);
    xq = linspace(pp.t(1),pp.t(end),100000);
    ctrlInt(1,:) = interp1(t,ctrl(1,:),xq,'previous');
    ctrlInt(2,:) = interp1(t,ctrl(2,:),xq,'next');
    ctrlInt(3,:) = interp1(t,ctrl(3,:),xq,'next');
    ctrlNorm = normOfVec(ctrlInt);
    plot(xq,ctrlInt','LineWidth',2)
    hold on
    plot(xq,ctrlNorm,'LineWidth',2,'color','k')
    ylabel('$a$ [mm/s2]')
end
if ~pp.cislunar
    xlabel('Number of orbits to TCA [-]')
else
    xlabel('Time to TCA [days]')
end
legend('R','T','N','$|\cdot|$','interpreter','latex')
grid on
axis tight
set(gca,'xdir','reverse')
hold off

% Ellipse B-plane
for k = 1:pp.n_conj
    xb     = xBall(:,k);
    x      = xMan(:,k);
    x_s    = x_sTCA(:,k);
    e2b    = eci2Bplane(xb(4:6),x_s(4:6));
    e2b    = e2b([1 3],:);
    PB     = e2b*P(:,:,k)*e2b';
    [semiaxes,cov2b] = defineEllipsoid(PB,smdLim(k));
    a          = semiaxes(1)*Lsc;
    b          = semiaxes(2)*Lsc;
    tt         = 0:0.001:2*pi;
    xx         = a*cos(tt);
    yy         = b*sin(tt);
    xx1        = mdLim*Lsc*cos(tt);
    yy1        = mdLim*Lsc*sin(tt);
    ellCov     = [xx; yy];
    ellB       = nan(2,length(tt));
    for j = 1:length(tt)
        ellB(:,j) = cov2b*ellCov(:,j);
    end
    figure()
    hold on    
    pOldB = e2b*(xb(1:3)-x_s(1:3))*Lsc;
    pNewB = e2b*(x(1:3)-x_s(1:3))*Lsc;
    plot(ellB(1,:),ellB(2,:),'k');
    plot(pNewB(1,:),pNewB(2,:),'o','LineWidth',2);
    plot(pOldB(1,:),pOldB(2,:),'k','marker','diamond');
    grid on 
    xlabel('$\xi$ [km]')
    ylabel('$\zeta$ [km]')
    hold off
    axis equal
end

end