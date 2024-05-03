clear
close all
for j = 2:9
    load(['SimOutput\diffManTimes\diffManTimesord',num2str(j),'.mat'])
    poc(:,j-1) = PoC;
    dvNorm(:,j-1) = normOfVec(dvs);
    compTime(:,j-1) = simTime;
end


figure
semilogy(tMan,poc)
xlabel('Orbits to TCA [-]')
ylabel('PoC after maneuver [-]')
grid on
legend('$n=2$','$n=3$','$n=4$','$n=5$','$n=6$','$n=7$','$n=8$','$n=9$','Interpreter','Latex')

figure
semilogy(tMan,dvNorm*pp.Vsc*1e6)
xlabel('Orbits to TCA [-]')
ylabel('$||\Delta v||$ [mm/s]')
grid on
legend('$n=2$','$n=3$','$n=4$','$n=5$','$n=6$','$n=7$','$n=8$','$n=9$','Interpreter','Latex')

figure
semilogy(tMan,compTime,'.')
xlabel('Orbits to TCA [-]')
ylabel('Computation Time [-]')
grid on
legend('$n=2$','$n=3$','$n=4$','$n=5$','$n=6$','$n=7$','$n=8$','$n=9$','Interpreter','Latex')
