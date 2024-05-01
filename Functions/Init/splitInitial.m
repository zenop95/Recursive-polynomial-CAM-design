function pp = splitInitial(pp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

s = pp.secondary;
n = pp.gmmOrder;
%% Transform into ECI
for k = 1:length(s)
    if pp.singleObject && k > 1
        tca = s(k).tca;
        pp.secondary((k-1)*n+1:k*n) = pp.secondary((k-2)*n+1:(k-1)*n);
        for j = 1:pp.gmmOrder
            pp.secondary((k-1)*n+j).tca = tca;
            pp.secondary((k-1)*n+j).cdm = s(k).cdm;
        end
        continue
    end
    original = s(k);
    if isempty(s(k).x0) && ~isempty(s(k).relState)
        relState = s(k).relState;                                             % [-] (6,1) Relative secondary state in RTN at the time of conjunction
        stateNli = pp.primary.cart0; %%%%%%%%%%%%%%%%%%%%%%%%% da sistemare perché va fatto con il secondario
        [r2ep,wp] = rtn2eci(pp.primary.cart0(1:3,1),pp.primary.cart0(4:6,1));             % [-] (3,3) DCM RTN of primary to ECI for the time of conjunction
        R2Ep      = [r2ep zeros(3); zeros(3) r2ep];
        state     = stateNli + R2Ep*relState;
        [r2es,ws] = rtn2eci(state(1:3),state(4:6));                        % [-] (3,3) DCM RTN of secondary to ECI for the time of conjunction
        R2Es      = [r2es zeros(3); zeros(3) r2es];
    elseif ~isempty(s(k).x0) && isempty(s(k).relState)
        state = s(k).x0;                                                   % [-] (6,1) Absolute secondary state in ECI at the time of conjunction
        [r2es,ws] = rtn2eci(state(1:3),state(4:6));                        % [-] (3,3) DCM RTN of secondary to ECI for the time of conjunction
        R2Es      = [r2es zeros(3); zeros(3) r2es];
        
        stateNli  = state;
    else
        error('either the relative or the absolute state of the secondary must be defined');
    end
    C0  = R2Es*s(k).C0*R2Es';

    %% Find direction of maximum nonlinearity
    fid = fopen('initial_state.dat', 'w');
    fprintf(fid, '%40.12f\n', pp.N*pp.dt);
    fprintf(fid, '%40.12f\n', pp.et);
    for i = 1:6 
        fprintf(fid, '%40.12f\n', stateNli(i));
    end
    fclose(fid);
    !wsl ./CppExec/nliGmm
    nli = load('trustRegion.dat');

    %% Perform the split along the found direction
    dirRtn = normalize(nli.*diag(s(k).C0),'norm');
    dir    = R2Es*dirRtn;
    xc     = state;
    C0c    = C0;
    [xi, Ci, wi] = VittaldevSplit(xc,C0c,dir,n);
    for j = 1:n
        secInd                  = j + (k-1)*n;           
        pp.secondary(secInd)    = original;
        pp.secondary(secInd).w  = wi(j);
        pp.secondary(secInd).x0 = xi(:,j);                     % state in ECI
        if isempty(s(k).x0) && ~isempty(s(k).relState)
            pp.secondary(secInd).relState = R2Ep'*(pp.secondary(secInd).x0 - pp.primary.cart0);
            pp.secondary(secInd).x0 = [];                     % state in ECI
        end
        pp.secondary(secInd).C0           = R2Es'*Ci(:,:,j)*R2Es;        % covariance in RTN of secondary
    end
end
end