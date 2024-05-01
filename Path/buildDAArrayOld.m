function J = buildDAArray(C,E,vars,order)

% vars to exclude
if     vars(1) == 1 && vars(end) == 3; varsNo = 4:6; 
elseif vars(1) == 4;                   varsNo = 1:3;
else;                                  varsNo = []; 
end

C(sum(E(:,varsNo),2)>0) = [];         
E(sum(E(:,varsNo),2)>0,:) = [];     
C(sum(E(:,vars),2)~=order) = []; 
E(sum(E(:,vars),2)~=order,:) = [];
E = E(:,vars);
n = length(vars);
if order == 1; J = zeros(1,3); else; J = zeros(n*ones(1,order)); end

for i = 1:length(C)
    cond1  = any(ismember(E(i,:),order));
    cond2  = any(ismember(E(i,:),order-1))                            && cond1  == 0;
    cond3  = any(ismember(E(i,:),order-2)) && any(ismember(E(i,:),2)) && cond2  == 0 && cond1 == 0;
    cond4  = any(ismember(E(i,:),order-2)) && any(ismember(E(i,:),1)) && cond3  == 0 && cond2 == 0 && cond1 == 0;
    cond5  = any(ismember(E(i,:),order-3)) && any(ismember(E(i,:),3)) && ~any(ismember(E(i,:),2)) && cond4  == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond6  = any(ismember(E(i,:),order-3)) && any(ismember(E(i,:),2)) && cond5  == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond7  = any(ismember(E(i,:),order-4)) && any(ismember(E(i,:),4)) && ~any(ismember(E(i,:),3)) && ~any(ismember(E(i,:),2)) && cond6  == 0 && cond5 == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond8  = any(ismember(E(i,:),order-4)) && any(ismember(E(i,:),3)) && ~any(ismember(E(i,:),2)) && cond7  == 0 && cond6 == 0 && cond5 == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond9  = any(ismember(E(i,:),order-4)) && any(ismember(E(i,:),2)) && cond8  == 0 && cond7 == 0 && cond6 == 0 && cond5 == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond10 = any(ismember(E(i,:),order-5)) && any(ismember(E(i,:),5)) && cond9  == 0 && cond8 == 0 && cond7 == 0 && cond6 == 0 && cond5 == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond11 = any(ismember(E(i,:),order-5)) && any(ismember(E(i,:),4)) && ~any(ismember(E(i,:),2)) && cond10 == 0 && cond9 == 0 && cond8 == 0 && cond7 == 0 && cond6 == 0 && cond5 == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond12 = any(ismember(E(i,:),order-5)) && any(ismember(E(i,:),3)) && cond11 == 0 && cond10 == 0 && cond9 == 0 && cond8 == 0 && cond7 == 0 && cond6 == 0 && cond5 == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    cond13 = any(ismember(E(i,:),order-6)) && any(ismember(E(i,:),3)) && cond12 == 0 && cond11 == 0 && cond10 == 0 && cond9 == 0 && cond8 == 0 && cond7 == 0 && cond6 == 0 && cond5 == 0 && cond4 == 0 && cond3 == 0 && cond2 == 0 && cond1 == 0;
    conds  = [cond1,cond2,cond3,cond4,cond5,cond6,cond7,cond8,cond9,cond10,cond11,cond12,cond13];
    switch find(conds)
        case 1
            ord1 = order;   % highest order
            ord2 = 0;       % medium order
            ord3 = 0;       % lowest order
        case 2
            ord1 = order-1; % highest order
            ord2 = 1;       % medium order
            ord3 = 0;       % lowest order
        case 3
            ord1 = order-2; % highest order
            ord2 = 2;       % medium order
            ord3 = 1;       % lowest order
        case 4
            ord1 = order-2; % highest order
            ord2 = 1;       % medium order
            ord3 = 1;       % lowest order
        case 5
            ord1 = order-3; % highest order
            ord2 = 3;       % medium order
            ord3 = 0;       % lowest order
        case 6
            ord1 = order-3; % highest order
            ord2 = 2;       % medium order
            ord3 = 1;       % lowest order
        case 7
            ord1 = order-4; % highest order
            ord2 = 4;       % medium order
            ord3 = 0;       % lowest order
        case 8
            ord1 = order-4; % highest order
            ord2 = 3;       % medium order
            ord3 = 1;       % lowest order
        case 9
            ord1 = order-4; % highest order
            ord2 = 2;       % medium order
            ord3 = 0;       % lowest order
        case 10
            ord1 = order-5; % highest order
            ord2 = 5;       % medium order
            ord3 = 0;       % lowest order
        case 11
            ord1 = order-5; % highest order
            ord2 = 4;       % medium order
            ord3 = 1;       % lowest order
        case 12
            ord1 = order-5; % highest order
            ord2 = 3;       % medium order
            ord3 = 2;       % lowest order
        case 13
            ord1 = order-6; % highest order
            ord2 = 3;       % medium order
            ord3 = 3;       % lowest order
        otherwise
            error('Something wrong')
    end
    ord2(ord2 == ord1) = 0;
    ord3(ord3 == ord2) = 0;
    ord3(ord3 == ord1) = 0;
    ind1 = find(E(i,:) == ord1);
    ind2 = find(E(i,:) == ord2);
    ind3 = find(E(i,:) == ord3);
    perm = unique(perms([repmat(ind1,[1,ord1]),repmat(ind2,[1,ord2]),repmat(ind3,[1,ord3])]),'rows');
    k = size(perm,1);

    indices = cell(1, order);
    % Convert the permutation indices to linear indices using sub2ind
    for j = 1:order
        indices{j} = perm(:, j);
    end
    linear_indices = sub2ind(size(J), indices{:});
    
    % Assign the value to the corresponding elements in J
    J(linear_indices) = C(i)/k;
end

end