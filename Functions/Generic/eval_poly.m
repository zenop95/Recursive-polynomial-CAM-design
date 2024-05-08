function p = eval_poly(C,E,PoI,orderMax)
% PoI vettore riga     
C(sum(E,2)>orderMax)   = []; 
E(sum(E,2)>orderMax,:) = [];
p = [];
for i = 1:size(PoI,1)
    
    if isempty(C)
        p = [p;0];
    else
        PoImx = ones(length(C),1)*PoI(i,:);
        p = [ p; sum(C.*prod((PoImx.^E),2)) ];
    end
    
end