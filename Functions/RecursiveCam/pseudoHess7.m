function J = pseudoHess7(dv0, K)
n = length(dv0);
L = nan(n,n,n,n,n,n);
for j = 1:n
    for k = 1:n
        for u = 1:n
            for l = 1:n
                for f = 1:n
                    L(:,u,k,j,l,f) = squeeze(dv0'*K(:,:,u,k,j,l,f));
                end
            end
        end
    end
end
J = pseudoHess6(dv0, L);
end