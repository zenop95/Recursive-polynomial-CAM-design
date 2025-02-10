name = ["md05.mat","md15.mat","md25.mat","md35.mat"];
% name = ["pc4.mat","pc5.mat","pc6.mat","pc7.mat"];
percMax = 99;
res = struct();
for j = 1:4
    load(name(j))
    res.type           = 'MD';
    res.sims(j).lim     = 0.5+(j-1);
    % res.type           = 'PoC';
    % res.sims(j).lim    = 10^(-(4+j-1));
    % c                   = c(c>prctile(c,1));
    c = fracOrb';
    res.sims(j).tFire   = [c nan(1,2170-length(c))];
    c                   = lowerThanPrc(compTime,percMax)';
    res.sims(j).simTime = [c nan(1,2170-length(c))];
    c                   = lowerThanPrc(md,percMax)';
    res.sims(j).md      = [c nan(1,2170-length(c))];
    c                   = lowerThanPrc(poc,percMax)';
    res.sims(j).poc     = [c nan(1,2170-length(c))];
    dv                  = normOfVec(dvs);
    c                   = lowerThanPrc(dv,percMax);
    res.sims(j).dv      = [c nan(1,2170-length(c))];
end
clearvars -except res
save('resMdPp')
% save('resPocPp')



 % for k = 1:2170
    %     if data.("d^* [km]")(k) > 0.5+(j-1)
    %         compTime(k) = nan;
    %         dvs(:,k) = nan;
    %         fracOrb(k) = nan;
    %         md(k) = nan;
    %         poc(k) = nan;
    %         t(k) = nan;
    %         xF(:,k) = nan;
    %     end
    % end