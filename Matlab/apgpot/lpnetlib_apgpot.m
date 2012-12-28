clear all;
clc;
%December 13 2012
%Sets up an LP from LPNetlib and solves it with lpapgpot

%Load the indices of the standard form 
%problems in lpnetlib
load '../standard_form_indices';
problem_index = 5;
problem_uf_ix = st_ix(problem_index);
P = UFget(problem_uf_ix);

prob_name = [P.name, '_run'];
A = P.A;
%Problem parameters
m = size(A,1);
n = size(A,2);

b = P.b;
c = P.aux.c;
fprintf('Loaded problem %s, m:%i, n:%i\n',P.name,m,n);

%End of problem definition

pars.rho = n+sqrt(n);
% call to lpapgpot (which calls apgpot):
[x,f,R] = lpapgpot(A,b,c,pars);

% compare with result from linprog:
[x1,f1] = linprog(c,[],[],A,b,zeros(n,1),[]);

fprintf('Results (low precision): \n');
fprintf('  %-25s%7.3e \n',' APGPOT Optimal value:',f);
fprintf('  %-25s%7.3e \n','LINPROG Optimal value:',f1);

% run again, but with different parameters:
%pars.echo    = 0;     %no printing
%pars.beta    = 0.9;   %different step-size reduction
%pars.rho     = 9.6*n; %different potential parameter (arbitrarily chosen)
%pars.tol     = 1e-5;  %higher precision
%[x,f,R]      = lpapgpot(A,b,c,pars);
%
%fprintf('Results (high precision): \n');
%fprintf('  %-25s%7.3e \n',' APGPOT Optimal value:',f);
%fprintf('  %-25s%7.3e \n','LINPROG Optimal value:',f1);

