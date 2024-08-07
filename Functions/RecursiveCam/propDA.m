function [lim,coeff,timeSubtr,xTca,xRet0] = ...
              propDA(DAorder,u,scale,validateFlag,pp)
% propDA performs the DA propagation to build the NLP
% 
% INPUT:
%        pp = [struct] optimization paramters structure
% 
% OUTPUT:
%        lim       = [-] PoC limit
%        coeffPoC  = [struct] structure with polynomial coefficient
%        timeSubtr = [struct] time to subtract from execution time to
%                             disregard writing and reading times
%        xTca      = [struct] states at TCA
%        PoC0      = [struct] PoC of the ballistic trajectory
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------
n_conj     = pp.n_conj;
n_man      = pp.n_man;
pocType    = pp.pocType;
t          = pp.t;
et         = pp.et;
Lsc        = pp.Lsc;
mu         = pp.mu;
HBR        = pp.HBR;
x_pTCA     = pp.x_pTCA;
x_sTCA     = pp.x_sTCA;
x_ref      = pp.xReference;
P          = pp.P;
PoCLim     = pp.PoCLim;
N = length(pp.t);
u = reshape(u,[],1);
scale = reshape(scale,[],1);
xTca   = nan(6,pp.n_conj);
bb = tic;
% Write file to pass to C++
fid = fopen('write_read/initial_state.dat', 'w');
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
fprintf(fid, '%40.16f\n', pp.gravOrd);
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
    for j = 1:6
        fprintf(fid, '%40.16f\n', x_ref(j));
    end
    for k = 1:n_conj
        for j = 1:3 
            for i = 1:3 
                fprintf(fid, '%40.16f\n', P(j,i,k));
            end
        end
    end
    for i = 1:n_man 
        fprintf(fid, '%40.16f\n', pp.thrustMagnitude);
    end
    for i = 1:n_man 
        for j = 1:3 
            fprintf(fid, '%40.16f\n', pp.thrustDirections(j,i));
        end
    end
    fprintf(fid, '%2i\n', pp.flagCA);
    fprintf(fid, '%2i\n', pp.flagTanSep);
    fprintf(fid, '%2i\n', pp.flagAlt);
    fprintf(fid, '%2i\n', pp.flagReturn);
end
for i = 1:length(u) 
    fprintf(fid, '%40.16f\n', u(i));
end
for i = 1:N 
    fprintf(fid, '%40.16f\n', t(i));
    fprintf(fid, '%2i\n', pp.canFire(i));
    fprintf(fid, '%2i\n', pp.isConj(i));
    fprintf(fid, '%2i\n', pp.isRet(i));
end
fclose(fid);
timeSubtr1 = toc(bb);                                                           % Exlude writing time from computation time measure
xRet0 = [];
%% Run the C++ Executable to perform the DA propagation
if ~validateFlag
    !wsl ./CppExec/polyProp
elseif validateFlag
    !wsl ./CppExec/validatePoly
    lim=[];coeff=[];timeSubtr=[];
    x    = reshape(load("write_read/constPart.dat"),6,pp.n_conj+1);             % If validating we only care about the TCA positions
    if any(pp.isRet); xRet0 = x(:,end); end
    xTca  = x(:,1:end-1);
    return;
end
b = tic;

%% Extract output from propagation
a         = load("write_read/constPart.dat");                                                         
for k = 1:n_conj
    xTca(:,k) = a(1+(k-1)*6:6*k);                                               % [-] (6,n_conj) Constant part of the propagated state and control
end


coeff  = struct();
if ~validateFlag
    coeff  = LoadCOSY('write_read/constraints.dat', ...
                   (3-2*pp.fixedDir-pp.fixedMag)*pp.n_man,pp.n_constr,0);
end
timeSubtr = toc(b) + timeSubtr1 + load("write_read/timeOut.dat")/1000 ;         % Exclude reading time from computation time measure
lim = log10(PoCLim);
end