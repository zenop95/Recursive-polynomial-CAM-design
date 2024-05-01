function J = pseudoHess9(dv0, Q)
n = length(dv0);
R = nan(n,n,n,n,n,n,n,n);
for j = 1:n
    for k = 1:n
        for u = 1:n
            for l = 1:n
                for f = 1:n
                    for r = 1:n
                        for q = 1:n
                            R(:,u,k,j,l,f,r,q) = squeeze(dv0'*Q(:,:,u,k,j,l,f,r,q));
                        end
                    end
                end
            end
        end
    end
end
J = pseudoHess8(dv0, R);
end