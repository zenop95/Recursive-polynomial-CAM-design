
function DrawCube2(R,color)

% R matrix of interval bounds [3x2]

a =  [R(1,1);R(2,1);R(3,1)];
b =  [R(1,2);R(2,2);R(3,1)];
c =  [R(1,2);R(2,2);R(3,2)];

plot3([a(1) b(1) b(1) a(1) a(1)],[a(2) a(2) b(2) b(2) a(2)],...
      [a(3) a(3) a(3) a(3) a(3)],color,'linewidth', 1);
hold on

plot3([a(1) b(1) b(1) a(1) a(1)],[a(2) a(2) b(2) b(2) a(2)],...
    [c(3) c(3) c(3) c(3) c(3)],color,'linewidth', 1);

plot3([a(1) a(1)],[a(2) a(2)],[a(3) c(3)],color,'linewidth', 1);
plot3([b(1) b(1)],[a(2) a(2)],[a(3) c(3)],color,'linewidth', 1);
plot3([b(1) b(1)],[b(2) b(2)],[b(3) c(3)],color,'linewidth', 1);
plot3([a(1) a(1)],[b(2) b(2)],[a(3) c(3)],color,'linewidth', 1);