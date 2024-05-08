syms x1 x2 x3 real;
n = 1e1;

p = SymDAPoly(coeff,[x1 x2 x3],[4 5 6],100);
span = linspace(-1e-2,1e-2,n);
val = nan(n,n);
x3   = 0;
for i = 1:n
    x1 = span(i);
    for j = 1:n
        x2 = span(j);
        val(i,j) = subs(p);
    end
end
[X,Y] = meshgrid(span,span);
surf(X,Y,val)
xlabel('R')
ylabel('T')
zlabel('PoC')
shading interp
colorbar