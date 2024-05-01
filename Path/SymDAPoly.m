function poly = SymDAPoly(coeff,x,vars,scaleG)

% eliminate DA variable that are not used (position)
C = coeff.C;
E = coeff.E;
C(sum(E(:,1:3),2)>0) = [];         
E(sum(E(:,1:3),2)>0,:) = [];     
E = E(:,vars);
c = nan(length(C),1);

syms polyn;
poly = 0*x(1);
for i = 1:length(C)
    c(i) = C(i)/scaleG^sum(E(i,:));
    poly = poly + c(i)*prod(x.^E(i,:));
end


end