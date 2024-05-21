function [] = plotOrbitCislunar(pp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n = 5000;
t = linspace(0,2*pi,n);
xx(:,1) = [pp.x_pTCA(1:2); pp.x_pTCA(4:5)];
xn(:,1) = [pp.x_pTCA(1:2); pp.x_pTCA(4:5)];

% for i = 2:n
%     xx(:,i) = propCr3bp(xx(:,i-1),zeros(3,1),t(i)-t(i-1),pp.mu);
% end
for i = 2:n
    xx(:,i) = propCr3bp2dNoCtrl(xx(:,i-1),t(i)-t(i-1),pp.mu);
end
for i = 2:n
    xn(:,i) = propCr3bp2dNoCtrl(xx(:,i-1),t(i)-t(i-1),pp.mu);
end
figure
plot(xx(1,:),xx(2,:),'-')
hold on
plot(-pp.mu,0,'o','linewidth',2)
plot(1-pp.mu,0,'o','linewidth',2)
hold off
xlabel('X [-]')
ylabel('Y [-]')
zlabel('Z [-]')
grid on
box on
end