function J = pseudoHess5(dv0,H)
n = length(dv0);
F = nan(n,n,n,n);
for j = 1:n
    for k = 1:n
        for u = 1:n
            F(:,u,k,j) = squeeze(dv0'*H(:,:,u,k,j));
        end
    end
end
J = pseudoHess4(dv0, F);
end