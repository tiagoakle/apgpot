clear all;
addpath '..' %Make sure that the path of prox_pot_toy is available
%clc;

%December 15 2012
%Sets up a random LP to compare with prox_pot_toy

seed = 0;
randn('seed',seed);
rand('seed',seed);
m = 5;
n = 15;
A = randn(m,n);

x0    = rand(n,1);
b     = A*x0; %Find a rhs which makes the problem feasible

y0    = randn(m,1);
z0    = rand(n,1);
c     = A'*y0 + z0; %Generate a dual feasible

%XXX: When comparing to prox pot toy these parameters must be compatible
pars.rho = n+sqrt(n);
pars.accel = 1;
pars.maxit = 10000;
pars.tolnorm = 2;

% call to lpapgpot (which calls apgpot):
[x,f,R] = lpapgpot(A,b,c,pars);

% compare with result from linprog:
[x1,f1] = linprog(c,[],[],A,b,zeros(n,1),[]);

fprintf('Results (low precision): \n');
fprintf('  %-25s%7.3e \n',' APGPOT Optimal value:',f);
fprintf('  %-25s%7.3e \n','LINPROG Optimal value:',f1);

%Run prox_pot_toy to compare
prox_pot_toy

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


