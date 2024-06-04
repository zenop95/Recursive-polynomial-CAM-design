function [CGT,timeSubtr] = Cauchy_Green_prop(DAorder,u,pp)
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
fprintf(fid, '%40.16f\n', Lsc);
fprintf(fid, '%40.16f\n', mu);
for j = 1:6 
    fprintf(fid, '%40.16f\n', x_pTCA(j));
end
for i = 1:length(u) 
    fprintf(fid, '%40.16f\n', u(i));
end
for i = 1:N 
    fprintf(fid, '%40.16f\n', t(i));
    fprintf(fid, '%2i\n', pp.canFire(i));
    fprintf(fid, '%2i\n', pp.isRet(i));
end
fclose(fid);
aidaInit(pp);                                                                   % Initialize AiDA dynamics (not needeed in paper)
timeSubtr1 = toc(bb);                                                           % Exlude writing time from computation time measure
%% Run the C++ Executable to perform the DA propagation
!wsl ./CppExec/CGTProp
b = tic;

%% Extract output from propagation
CGT  = load("write_read/constPart.dat");                                                         

timeSubtr = toc(b) + timeSubtr1 + load("write_read/timeOut.dat")/1000 ;         % Exclude reading time from computation time measure

end