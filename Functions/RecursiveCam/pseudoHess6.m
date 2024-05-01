function J = pseudoHess6(dv0, L)
n = length(dv0);
H = nan(n,n,n,n,n);
for j = 1:n
    for k = 1:n
        for u = 1:n
            for l = 1:n
                H(:,u,k,j,l) = squeeze(dv0'*L(:,:,u,k,j,l));
            end
        end
    end
end
J = pseudoHess5(dv0, H);
end