placeFigure()
js = reshape(reshape(1:9,3,3)',1,[]);
for j = 1:9
    subplot(3,3,j)
    plot(1,Ys(js(j),1)*pp.Asc*pp.ctrlMax*1e6,'o') 
    hold on
    plot(2:11,Ys(js(j),2:11)*pp.Asc*pp.ctrlMax*1e6,'o') 
    plot(12:14,Ys(js(j),12:14)*pp.Asc*pp.ctrlMax*1e6,'o')
    plot(15:18,Ys(js(j),15:end)*pp.Asc*pp.ctrlMax*1e6,'o')
    grid on
end
% ylabel('$u$ [mm/s2]')
% xlabel('Iteration [-]')

figure()
    semilogy(1:10,er(1:10)*pp.Asc*pp.ctrlMax*1e6,'o') 
grid on
hold on
    semilogy(11:13,er(11:13)*pp.Asc*pp.ctrlMax*1e6,'o') 
    semilogy(14:17,er(14:17)*pp.Asc*pp.ctrlMax*1e6,'o') 
xlabel('Iteration [-]')
ylabel('$e$ [mm/s2]')
plot([1,sum(iters)],pp.tol*pp.Asc*pp.ctrlMax*1e6*ones(2,1),'k--')
hold off
axis tight

% figure()
% subplot(3,1,1)
% plot(Ys(1,:)*pp.Vsc*pp.ctrlMax*1e6,'.')
% grid on
% subplot(3,1,2)
% plot(Ys(2,:)*pp.Vsc*pp.ctrlMax*1e6,'.')
% grid on
% subplot(3,1,3)
% plot(Ys(3,:)*pp.Vsc*pp.ctrlMax*1e6,'.')
% grid on
% grid on
% ylabel('$\Delta v$ [mm/s]')
% xlabel('Iteration [-]')
% 
% figure()
% semilogy(reshape(er,1,[])*pp.Vsc*pp.ctrlMax*1e6,'.')
% grid on
% hold on
% xlabel('Iteration [-]')
% ylabel('$e$ [mm/s]')
% plot([1,sum(iters)],pp.tol*pp.Vsc*pp.ctrlMax*1e6*ones(2,1),'k--')
% yticks([1e-1,1e1, 1e3])
% hold off
% axis tight
