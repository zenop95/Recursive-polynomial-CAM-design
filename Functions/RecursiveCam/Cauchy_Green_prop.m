function [maxEig,timeSubtr] = Cauchy_Green_prop(DAorder,u,pp)
% Cauchy_Green_prop performs the DA propagation
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
n_man      = pp.n_man;
t          = pp.t;
et         = pp.et;
Lsc        = pp.Lsc;
mu         = pp.mu;
N = length(pp.t);
u = reshape(u,[],1);
bb = tic;
% Write file to pass to C++
fid = fopen('write_read/initial_state.dat', 'w');
fprintf(fid, '%2i\n',     N);
fprintf(fid, '%2i\n',     n_man);
fprintf(fid, '%2i\n',     pp.m);
fprintf(fid, '%2i\n',     pp.lowThrust);
fprintf(fid, '%2i\n',     DAorder);
fprintf(fid, '%40.16f\n', et);
fprintf(fid, '%40.16f\n', pp.primary.T);
fprintf(fid, '%40.16f\n', Lsc);
fprintf(fid, '%40.16f\n', mu);
for j = 1:6 
    fprintf(fid, '%40.16f\n', pp.x_pTCA(j));
end
for i = 1:N 
    fprintf(fid, '%40.16f\n', t(i));
end
fclose(fid);
timeSubtr1 = toc(bb);                                                           % Exlude writing time from computation time measure
%% Run the C++ Executable to perform the DA propagation
!wsl ./CppExec/CGTProp
b = tic;

%% Extract output from propagation
s    = load("write_read/constPart.dat");
STMs  = reshape(s,6,6,N-1);
for i = 1:N-1
    STM = STMs(:,:,i);
    CGT = STM'*STM; CGT = CGT(4:6,4:6);
    [evects,D] = eig(CGT);
    evals(:,i) = diag(D);
    [c(i),a(i)] = min(evals(:,i));
    maxEig(:,i) = evects(:,a(i));
end
timeSubtr = toc(b) + timeSubtr1 + load("write_read/timeOut.dat")/1000 ;         % Exclude reading time from computation time measure

end