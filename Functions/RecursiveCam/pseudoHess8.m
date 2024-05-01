function J = pseudoHess8(dv0, R)
n = length(dv0);
K = nan(n,n,n,n,n,n,n);
for j = 1:n
    for k = 1:n
        for u = 1:n
            for l = 1:n
                for f = 1:n
                    for r = 1:n
                        K(:,u,k,j,l,f,r) = squeeze(dv0'*R(:,:,u,k,j,l,f,r));
                    end
                end
            end
        end
    end
end
J = pseudoHess7(dv0, K);
end