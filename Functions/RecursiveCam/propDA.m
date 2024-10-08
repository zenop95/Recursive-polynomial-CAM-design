function [lim,coeff,timeSubtr,xTca,xRet0,x_sec,deltaTca] = ...
              propDA(DAorder,u,validateFlag,pp)
% propDA performs the DA propagation to build the NLP
% 
% INPUT:
%        DAorder      = [-]    (1,1) Maximum order of the polynomial constraints
%        u            = [-]    (m,1) Nominal control (usually 0)
%        validateFlag = [bool] (1,1) Toggle validation propagation
%        pp           = [struct] optimization paramters structure
% 
% OUTPUT:
%        lim       = [-]      (1,1)      PoC limit
%        coeff     = [struct]            Structure with polynomial coefficients
%        timeSubtr = [s]      (1,1)      Time to subtract from execution time to
%                                        disregard writing and reading times
%        xTca      = [-]      (6,n_conj) States of the primary at TCAs 
%        xRet0     = [-]      (6,1)      State at SK node
%        x_sec     = [-]      (6,n_conj) States of the secondary at TCAs 
%        deltaTca  = [-]      (n_conj,1) TCA shift due to maneuver
%
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------------
%% initialize variables from input structure
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
PoCLim     = pp.PoCLim;
N          = length(pp.t);
u          = reshape(u,[],1);
xTca       = nan(6,pp.n_conj);
bb         = tic;
%% Write file to pass to C++
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
fprintf(fid, '%2i\n',     pp.gravOrd);
% fprintf(fid, '%2i\n',     2);
fprintf(fid, '%40.16f\n', pp.ctrlMax);
fprintf(fid, '%2i\n',     pp.flagMd);
fprintf(fid, '%2i\n',     pp.flagPoCTot);
fprintf(fid, '%40.16f\n', pp.primary.n);
for k = 1:n_conj
    fprintf(fid, '%40.16f\n', pp.secondary(k).n);
end
for j = 1:6 
    fprintf(fid, '%40.16f\n', x_pTCA(j));
end
for k = 1:n_conj
    for j = 1:6
        fprintf(fid, '%40.16f\n', x_sTCA(j,k));
    end
end
if ~validateFlag
    for k = 1:n_conj
        fprintf(fid, '%40.16f\n', HBR(k));
    end
    for k = 1:n_conj
        for j = 1:6 
            for i = 1:6 
                fprintf(fid, '%40.16f\n', pp.Cp(j,i,k));
            end
        end
        for j = 1:6 
            for i = 1:6 
                fprintf(fid, '%40.16f\n', pp.Cs(j,i,k));
            end
        end
    end
    for i = 1:n_man 
        for j = 1:3 
            fprintf(fid, '%40.16f\n', pp.thrustDirections(j,i));
        end
    end
    fprintf(fid, '%2i\n', pp.flagCA);
    fprintf(fid, '%2i\n', pp.flagReturn);
    fprintf(fid, '%2i\n', pp.flagErrReturn);
    fprintf(fid, '%2i\n', pp.flagMeanSma);
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
    out = reshape(load("write_read/constPart.dat"),6,2*pp.n_conj+1);            
    xTca = out(:,1:n_conj);
    x_sec = out(:,n_conj+1:2*n_conj);
    deltaTca = load("write_read/tcaOut.dat")*pp.Tsc;             
    if any(pp.isRet); xRet0 = out(:,end); end
    return;
end
b = tic;

%% Extract output from propagation
a = load("write_read/constPart.dat");                                                         
for k = 1:n_conj
    xTca(:,k) = a(1+(k-1)*6:6*k);                                               % [-] (6,n_conj) Constant part of the propagated state and control
end
if pp.flagReturn || pp.flagErrReturn || pp.flagMeanSma
    xRet0 = a(end-5:end);                                               % [-] (6,n_conj) Constant part of the propagated state and control
end

if pp.flagMd
    lim = pp.mdLim;
else
    lim = log10(PoCLim);
end


coeff  = struct();
if ~validateFlag
    coeff  = LoadCOSY('write_read/constraints.dat', ...
                   (3-2*pp.fixedDir)*pp.n_man,pp.n_constr,0);
end

timeSubtr = toc(b) + timeSubtr1 + load("write_read/timeOut.dat")/1000 ;         % Exclude reading time from computation time measure
end