function [lim,smdLim,coeffPoC,timeSubtr,xTca,metric,dist,convRadius] = ...
              propDAPoly(DAorder,u,scale,validateFlag,convRadFlag,pp)

n_conj     = pp.n_conj;
n_man      = pp.n_man;
pocType    = pp.pocType;
t          = pp.t;
et         = pp.et;
mdLim      = pp.mdLim;
Lsc        = pp.Lsc;
mu         = pp.mu;
HBR        = pp.HBR;
x_pTCA     = pp.x_pTCA;
x_sTCA     = pp.x_sTCA;
P          = pp.P;
PoCLim     = pp.PoCLim;
convRadius = nan;
N = length(pp.t);
u = reshape(u,[],1);
scale = reshape(scale,[],1);
xTca   = nan(6,pp.n_conj);
smdLim = nan(pp.n_conj,1);
bb = tic;
fid = fopen('initial_state.dat', 'w');
fprintf(fid, '%2i\n',     N);
fprintf(fid, '%2i\n',     n_conj);
fprintf(fid, '%2i\n',     n_man);
fprintf(fid, '%2i\n',     pp.m);
fprintf(fid, '%2i\n',     pp.cislunar);
fprintf(fid, '%2i\n',     pp.lowThrust);
fprintf(fid, '%2i\n',     DAorder);
fprintf(fid, '%2i\n',     pocType);
fprintf(fid, '%40.16f\n', et);
fprintf(fid, '%40.16f\n', Lsc);
fprintf(fid, '%40.16f\n', mu);
for j = 1:6 
    fprintf(fid, '%40.16f\n', x_pTCA(j));
end
if ~validateFlag
    for k = 1:n_conj
        fprintf(fid, '%40.16f\n', HBR(k));
    end
    for k = 1:n_conj
        for j = 1:6
            fprintf(fid, '%40.16f\n', x_sTCA(j,k));
        end
    end
    for k = 1:n_conj
        for j = 1:3 
            for i = 1:3 
                fprintf(fid, '%40.16f\n', P(j,i,k));
            end
        end
    end
    for i = 1:length(scale) 
        fprintf(fid, '%40.16f\n', scale(i));
    end
    for i = 1:n_man 
        fprintf(fid, '%40.16f\n', pp.thrustMagnitude);
    end
    for i = 1:n_man 
        for j = 1:3 
            fprintf(fid, '%40.16f\n', pp.thrustDirections(j,i));
        end
    end
end
for i = 1:length(u) 
    fprintf(fid, '%40.16f\n', u(i));
end
for i = 1:N 
    fprintf(fid, '%40.16f\n', t(i));
    fprintf(fid, '%2i\n', pp.canFire(i));
    fprintf(fid, '%2i\n', pp.isConj(i));
end
fclose(fid);
aidaInit(pp,'primary');
timeSubtr1 = toc(bb);
% Run the C++ Executable to perform the DA propagation
if ~validateFlag && ~convRadFlag
    !wsl ./CppExec/polyProp
elseif validateFlag && ~convRadFlag
    !wsl ./CppExec/validatePoly
elseif ~validateFlag && convRadFlag
    !wsl ./CppExec/polyConvRadius
    convRadius = load('convRad.dat');
    convRadius = reshape(convRadius,[],N);
else
    error('Invalid flag combination')
end
if convRadFlag; lim=[];smdLim=[];coeffPoC=[];timeSubtr=[];xTca=[];metric=[];dist=[]; 
    return; 
elseif validateFlag
    lim=[];smdLim=[];coeffPoC=[];timeSubtr=[];dist=[];metric=[];
    xTca = reshape(load("constPart.dat"),6,pp.n_conj);                          % [-] (6,1) Constant part of the propagated state and control
    return; 
end
b = tic;
% Polynomials extraction 
a         = load("constPart.dat");                                                         
for k = 1:n_conj
    xTca(:,k) = a(1+(k-1)*6:6*k);                                               % [-] (6,n_conj) Constant part of the propagated state and control
end
metric = a(n_conj*6 + 1);                                                       % [-] (1,1) Collision metric with no maneuver
detPB  = a(n_conj*6 + 1 + (1:k));                                               % [-] (1,1) Determinant of the combined covariance at TCA (2d)
P_B    = reshape(a(end-n_conj*4+1:end),2,2,n_conj);                             % [-] (1,1) Determinant of the combined covariance at TCA (2d)
timeSubtr = toc(b) + timeSubtr1 + load("timeOut.dat")/1000 ;
for k = 1:n_conj
    switch pocType
        case 0
            smdLim(k)   = -2*log(2*PoCLim*sqrt(detPB(k))/HBR(k)^2);                                  % [-] (1,1) SMD limit computed with Alfriend and Akella's formula applied to PC
        case 1
            smdLim(k)   = PoC2SMD(P_B(:,:,k), HBR(k), PoCLim, 3, 1, 1e-3, 200);                          % [-] (1,1) SMD limit computed with Chan's formula applied to PC
        otherwise
            error('invalid PoC type')
    end
end
lim = log10(PoCLim);
coeffPoC  = struct();
if ~validateFlag && ~convRadFlag
    coeffPoC  = LoadCOSY('metricPoly.dat',(3-2*pp.fixedDir-pp.fixedMag)*pp.n_man,1,0);
elseif validateFlag
end
end