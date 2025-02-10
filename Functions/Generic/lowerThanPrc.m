function [out] = lowerThanPrc(data,prc)
    out = data(data<prctile(data,prc));
end