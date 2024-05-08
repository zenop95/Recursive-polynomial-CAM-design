function DAArray = buildDAArray(C,E,n)
% computeDvGreedy Solves the polynomial CAM optimization problem using a
% recursive greedy solution
%
% INPUT:  C = [-] Vector with values of the coefficients of the expansion
%         E = [-] Matrix with positions of the coefficients of the expansion
%         n = [-] Expansion order of the associated polynomial truncation
%
% OUTPUT: DAArray = [cell] Cell array for nth-order tensor
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------

C(sum(E,2)~=n)   = [];                                                          % [-] (*,*) Remove coefficients of different order than n
E(sum(E,2)~=n,:) = [];                                                          % [-] (*,*) Clean-up E matrix to be compatible with the new C
m = size(E,2);                                                                  % [-] (1,1) Number of coefficients
if n == 1; DAArray = zeros(1,m); else; DAArray = zeros(m*ones(1,n)); end        % [-] (1,m) Initialize gradient as row vector in the case of n = 1

% Each coefficient will appear more than once in the high-order tensor
% because of the super-symmetry. For example, in the Hessian matrix, the
% diagoanl coefficients only appear once, but the off-diagonal ones appear
% twice. The original coefficients C include the contribution of each occurrence,
% so they must be divided by the number of occurrences.
for i = 1:length(C)
    indR  = cell(m,1);                                                          % [cell] (m,1) Initialize indeces
    for j = 1:m
        indR{j} = repmat(j,1,E(i,j));                                           % [cell] (1,m) Build indeces of symmetric
    end
    perm = unique(perms(cat(2,indR{:})),'rows');                                % [-] (1,m) Compute all possible permutations of the indices to build super-symmetric tensor
    k = size(perm,1);                                                           % [-] (1,m) Number of possible permutations (number of occurrences of C in the high-order tensor)

    indices = cell(1, n);                                                       % Initialize indices cell array for the coefficient C
    for j = 1:n
        indices{j} = perm(:, j);                                                % All possible combinations of the indeces
    end
    linear_indices = sub2ind(size(DAArray), indices{:});                        % [-] Convert the tensor indices to linear indeces to associate with coefficients
    DAArray(linear_indices) = C(i)/k;                                           % [-] Assign value scaled by number of permutations because the coefficient C includes the contribution of every permutation
end

end