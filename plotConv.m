placeFigure()
js = reshape(reshape(1:9,3,3)',1,[]);
for j = 1:9
    subplot(3,3,j)
    plot(Ys(js(j),:)*pp.Asc*pp.ctrlMax*1e6,'.')
    grid on
    % ylabel('$\Delta v$ [mm/s]')
end
ylabel('$u$ [mm/s2]')
xlabel('Iteration [-]')

figure()
semilogy(reshape(er,1,[])*pp.Vsc*pp.ctrlMax*1e6,'.')
grid on
hold on
xlabel('Iteration [-]')
ylabel('$e$ [mm/s2]')
plot([1,sum(iters)],pp.tol*pp.Vsc*pp.ctrlMax*1e6*ones(2,1),'k--')
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
