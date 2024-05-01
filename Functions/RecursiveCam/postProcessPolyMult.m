function [] = postProcessPolyMult(smdLim,xBall,x,dv,distVal,metricVal,simTime,M,ns,pp)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

mdLim      = pp.mdLim;
Lsc        = pp.Lsc;
Vsc        = pp.Vsc;
x_sTCA     = pp.x_sTCA;
P          = pp.P;

if pp.lowThrust
    a  = dv;
    dv = nan(3,pp.M);
    for j = 1:pp.M
        dv(:,j) = squeeze(sum(a(:,:,j).*repmat([diff(pp.t),0],3,1),2));
    end
end

% Ellipse B-plane
colors = [0 0.4470 0.7410];
e2b = eci2Bplane(xBall(4:6,1),x_sTCA(4:6));
e2b = e2b([1 3],:);
PB  = e2b*P*e2b';
[semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
a          = semiaxes(1)*Lsc;
b          = semiaxes(2)*Lsc;
tt         = 0:0.001:2*pi;
xx         = a*cos(tt);
yy         = b*sin(tt);
xx1        = mdLim*Lsc*cos(tt);
yy1        = mdLim*Lsc*sin(tt);
ellCov     = [xx; yy];
ellB       = nan(2,length(tt));
for k = 1:length(tt)
    ellB(:,k) = cov2b*ellCov(:,k);
end
figure
hold on    
pOldB = e2b*(xBall(1:3)-x_sTCA(1:3))*Lsc;
for j = 1:M
    pNewB(:,j) = e2b*(x(1:3,j)-x_sTCA(1:3))*Lsc;
end
plot(ellB(1,:),ellB(2,:),'k');
% plot(xx1,yy1,'k--');
plot(pOldB(1),pOldB(2),'color',colors, ...
    'Marker','diamond','HandleVisibility','off')
s = scatter(pNewB(1,:),pNewB(2,:),[],'filled');
s.SizeData = 20;
colormap jet;
% plot(pNewB(1,:),pNewB(2,:),'o');
grid on 
xlabel('$\xi$ [km]')
ylabel('$\zeta$ [km]')
hold off
axis equal
cb = colorbar;
ylabel(cb,'Orbits to TCA [-]','Interpreter','latex')

% figure
% plot(1:ordMax,metricVal)
% grid on
% xlabel('DA order [-]')
% ylabel('PoC [-]')

% if length(ns) > 1
%     figure
%     plot(ns,simTime)
%     grid on
%     xlabel('Orbits to TCA [-]')
%     ylabel('$T_{sim}$ [s]')

%     figure
%     yyaxis left
%     semilogy(ns,norm(distVal(1:3,:))*pp.Lsc)
%     grid on
%     xlabel('Orbits to TCA [-]')
%     ylabel('Position target error [km]')
%     yyaxis right
%     semilogy(ns,norm(distVal(4:6,:))*pp.Vsc*1e3)
%     grid on
%     ylabel('Velocity target error [m/s]')

%     figure
%     semilogy(ns,metricVal)
%     grid on
%     xlabel('Orbits to TCA [-]')
%     ylabel('Validated $P_C$ [-]')
    % 
    % figure
    % plot(ns,smd)
    % hold on
    % plot([ns(1),ns(end)],smdLim*[1,1],'k--')
    % grid on
    % xlabel('Orbits to TCA [-]')
    % ylabel('$d_m^2$ [-]')
    % figure
    % plot(ns,sqrt(smd)*Lsc)
    % hold on
    % plot([ns(1),ns(end)],mdLim*Lsc*[1,1],'k--')
    % grid on
    % xlabel('Orbits to TCA [-]')
    % ylabel('$d_{miss}$ [-]')
    
%     figure
%     plot(ns,dv'*Vsc*1e6)
%     hold on
%     plot(ns,normOfVec(dv)*Vsc*1e6,'k--')
%     hold off
%     grid on
%     xlabel('Orbits to TCA [-]')
%     ylabel('$\Delta V$ [mm/s]')
%     legend('R','T','N','||\cdot||')
end