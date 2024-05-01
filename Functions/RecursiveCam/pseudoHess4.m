function J = pseudoHess4(dv0, F)
n = length(dv0);
for j = 1:n
    for k = 1:n
        G(:,k,j) = squeeze(dv0'*F(:,:,k,j));
    end
end
J = pseudoHess3(dv0, G);
end