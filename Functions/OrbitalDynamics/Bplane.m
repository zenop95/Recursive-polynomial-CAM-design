function [Pb,rb,smd] = Bplane(x_p,x_s,P)

e2b    = eci2Bplane(x_p(4:6),x_s(4:6));
e2b    = e2b([1 3],:);
Pb     = e2b*P*e2b';
rb     = e2b*(x_p(1:3)-x_s(1:3));
smd    = dot(rb,Pb\rb);

end