function [] = postProcess(xBall,xManTca,xManSec,xRetMan,xRetBall,lim,ctrl,deltaTca,simTime,pp)
% postProcess plots the relevant data 
% 
% INPUT: 
%        xBall
%        xMan
%        lim
%        ctrl
%        simTime
%        pp = [struct] optimization paramters structure
% 
% OUTPUT:
%        pp    = [struct] optimization paramters structure
%        t_man = [-] Maneuvering times in orbit units
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

Lsc        = pp.Lsc;
Vsc        = pp.Vsc;
P          = pp.P;
pp.t       = -pp.t;
%nodes for low-thrust
if pp.lowThrust
    t_lt            = pp.t; 
%     t_lt(pp.isConj) = [];
    dt_lt           = diff(t_lt); 
    dt_lt           = abs(dt_lt(pp.canFire));
    s               = strfind(pp.canFire',[1 0]) + 1; % Nodes in which the thruster is shut-off
    ttt = [];
    for j = 1:length(s)
        ttt = [ttt; pp.t(s(j))-1e-6];
    end
end
% compute PoC after maneuver
PoC = nan(pp.n_conj,1);
for k = 1:pp.n_conj
    x      = xManTca(:,k);
    x_s    = xManSec(:,k);
    xb     = xBall(:,k);
    STMp   = CWStateTransition(pp.primary.n^(3/2),deltaTca(k)/pp.Tsc,0,1);
    STMs   = CWStateTransition(pp.secondary(k).n^(3/2),deltaTca(k)/pp.Tsc,0,1);
    Cpprop = STMp*pp.Cp(:,:,k)*STMp';
    Csprop = STMs*pp.Cs(:,:,k)*STMs';
    Pp     = Cpprop(1:3,1:3);
    Ps     = Csprop(1:3,1:3);
    r2ep   = rtn2eci(x(1:3),x(4:6));
    r2es   = rtn2eci(x_s(1:3),x_s(4:6));
    P      = r2ep*Pp*r2ep' + r2es*Ps*r2es';
    [PB,p,smd] = Bplane(x,x_s,P);
    md = norm(p)*Lsc;
    switch pp.pocType
    case 0
        PoC(k) = constantPc(p,PB,pp.HBR(k));                                   % [-] (1,1) PoC computed with Chan's formula
    case 1
        PoC(k) = poc_Chan(pp.HBR(k),PB,smd,3);                                 % [-] (1,1) PoC computed with Chan's formula
    case 2
        PoC(k) = maximumPc(p,PB,pp.HBR(k));                                   % [-] (1,1) PoC computed with Chan's formula
    case 3
        PoC(k) = norm(p)*pp.scaling(1);                                   % [-] (1,1) PoC computed with Chan's formula
    otherwise
        error('invalid PoC type')
    end
    if ~pp.flagMd
        switch pp.pocType
        case 0
            smdLim   = -2*log(2*pp.PoCLim*sqrt(det(PB))/pp.HBR(k)^2);        % [-] (1,1) SMD limit computed with Alfriend and Akella's formula
        case 1
            smdLim   = PoC2SMD(PB, pp.HBR(k), pp.PoCLim, 3, 1, 1e-3, 200);   % [-] (1,1) SMD limit computed with Chan's formula
        case 2
            smdLim   = pp.HBR(k)^2/(exp(1)*sqrt(det(PB))*pp.PoCLim);                         % [-] (1,1) SMD limit computed with Maximum formula
        otherwise
            error('invalid PoC type')
        end
    else
        smdLim   = pp.mdLim;                                               % [-]   (1,1) PoC limit;     % [-] (1,1) SMD limit computed with Miss distance
        PB       = eye(2);    
    end
    [semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
    a          = semiaxes(1)*Lsc;
    b          = semiaxes(2)*Lsc;
    tt         = 0:0.001:2*pi;
    xx         = a*cos(tt);
    yy         = b*sin(tt);
    ellCov     = [xx; yy];
    ellB       = nan(2,length(tt));
    for j = 1:length(tt)
        ellB(:,j) = cov2b*ellCov(:,j);
    end
    figure('Renderer', 'painters', 'Position', [600*(k-1) 600 560 300])
    hold on    
    e2b    = eci2Bplane(x(4:6),x_s(4:6));
    e2b    = e2b([1 3],:);
    pOldB = e2b*(xb(1:3)-pp.x_sTCA(1:3,k))*Lsc;
    plot(ellB(2,:),ellB(1,:),'k');
    plot(p(2,:)*Lsc,p(1,:)*Lsc,'o','LineWidth',2);
    plot(pOldB(2,:),pOldB(1,:),'k','marker','diamond');
    grid on 
    xlabel('$\zeta$ [km]')
    ylabel('$\xi$ [km]')
    hold off
    axis equal
    box on
end
poc_tot = PoCTot(PoC);

% Validate return
if pp.flagReturn || pp.flagErrReturn || pp.flagMeanSma
    errRet = xRetMan - xRetBall;
    finalCoeMan   = cartesian2kepler(xRetMan,1);     
    finalCoeBall  = cartesian2kepler(xRetBall,1);     
    finalMeanCoeMan  = osculating2mean(finalCoeMan,1,pp.Lsc);
    finalMeanCoeBall = osculating2mean(finalCoeBall,1,pp.Lsc);
end

disp(['Solver: ', pp.solvingMethod])
disp(['Computation time ',num2str(simTime), ' s'])
n   = size(ctrl,2);
disp(['Number of conjunctions: ', num2str(pp.n_conj)])
disp(['tca shift = ', num2str(deltaTca(1)), ' s'])

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
disp(['PoC after validation ',num2str(poc_tot)]);
disp(['MD after validation ',num2str(md(1)), ' km']);
if pp.flagReturn || pp.flagErrReturn ||  pp.flagMeanSma
    disp(['Position error in return ',num2str(norm(errRet(1:3))*pp.Lsc*1e3), ' m']);
    disp(['Velocity error in return ',num2str(norm(errRet(4:6))*pp.Vsc*1e6), ' mm/s']);
    disp(['SMA error in return ',num2str((finalCoeMan.a-finalCoeBall.a)*pp.Lsc*1e3), ' m']);
    disp(['ECC error in return ',num2str(finalCoeMan.ecc-finalCoeBall.ecc)]);
    disp(['ARG error in return ',num2str(rad2deg(finalCoeMan.w-finalCoeBall.w)), ' deg']);
    disp(['RAAN error in return ',num2str(rad2deg(finalCoeMan.RAAN-finalCoeBall.RAAN)), ' deg']);
    disp(['TA error in return ',num2str(rad2deg(finalCoeMan.theta-finalCoeBall.theta)), ' deg']);
    disp(['INC error in return ',num2str(rad2deg(finalCoeMan.inc-finalCoeBall.inc)), ' deg']);
end
if pp.flagMeanSma
    disp(['Semi-major axis shift ',num2str((finalMeanCoeMan.a - finalMeanCoeBall.a)*pp.Lsc*1e3), ' m']);
    disp(['Eccentricity shift ',num2str((finalMeanCoeMan.ecc*finalMeanCoeMan.a - finalMeanCoeBall.ecc*finalMeanCoeBall.a)*pp.Lsc*1e3), ' m']);
end
if pp.flagMd
    disp(['Limit: ',num2str(sqrt(pp.mdLim)*pp.Lsc), ' km'])
else
    disp(['Limit: ',num2str(lim)])
end
ctrlNorm = normOfVec(ctrl);

figure('Renderer', 'painters', 'Position', [1000 300 560 150])
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
    % stem(t(ctrlNorm~=0),ctrlNorm(ctrlNorm~=0),'color','k','LineWidth',2)
    % plot(t(ctrlNorm==0),ctrlNorm(ctrlNorm==0),'color','k')
    ylabel('$\Delta v$ [mm/s]')
    % set(gca, 'XDir','reverse')
else
    ctrlN  = [];
    ctrlN1 = ctrl*pp.Asc*1e6;
    ccc = 0;
    for i = 1:pp.N
        if ~pp.canFire(i) 
            ctrlN = [ctrlN zeros(3,1)];
        else
            ccc = ccc + 1;
            ctrlN = [ctrlN ctrlN1(:,ccc)];
        end
    end
    [t,o] = sort([pp.t;ttt],'descend');
    if ~pp.cislunar
        t     = t/pp.T;
    else
        t  = t*pp.Tsc/86400;
    end
    ctrlN  = [ctrlN, ctrlN1];
    ctrlN  = ctrlN(:,o);
    xq = linspace(t(1),t(end),1000);
    ctrlInt(1,:) = interp1(t,ctrlN(1,:),xq,'previous');
    ctrlInt(2,:) = interp1(t,ctrlN(2,:),xq,'previous');
    ctrlInt(3,:) = interp1(t,ctrlN(3,:),xq,'previous');
    ctrlNorm = normOfVec(ctrlInt);
    plot(xq,ctrlInt','LineWidth',2)
    hold on
    plot(xq,ctrlNorm,'LineWidth',2,'color','k')
    ylabel('$a$ [mm/s2]')
end
if ~pp.cislunar
    xlabel('Number of orbits [-]')
else
    xlabel('Time to TCA [days]')
end
if pp.cislunar && ~pp.lowThrust; legend('$\Delta V_x$','$\Delta V_y$','$\Delta V_z$','$||\Delta V||$','interpreter','latex'); end
if pp.cislunar && pp.lowThrust; legend('u_x','u_y','u_z','$|\cdot|$','interpreter','latex'); end
grid on
axis tight
yl = ylim;
conjs = pp.t(pp.isConj)/pp.T;
for j = 1:pp.n_conj
    plot(conjs(j)*ones(1,2),yl,'k--')
end
legend('R','T','N','$|\cdot|$','interpreter','latex')
hold off
% saveas(gcf, 'dv', 'epsc') %Save figure

end