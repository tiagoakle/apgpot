% generate a random LP:
m  = 5;
n  = 15;
A  = randn(m,n);
xx = rand(n,1);
b  = A*xx;
ss = rand(n,1);
c  = A'*randn(m,1) + ss;

% call to lpapgpot (which calls apgpot):
[x,f,R] = lpapgpot(A,b,c);

% compare with result from linprog:
[x1,f1] = linprog(c,[],[],A,b,zeros(n,1),[]);

fprintf('Results (low precision): \n');
fprintf('  %-25s%7.3e \n',' APGPOT Optimal value:',f);
fprintf('  %-25s%7.3e \n','LINPROG Optimal value:',f1);

% run again, but with different parameters:
pars.echo    = 0;     %no printing
pars.beta    = 0.9;   %different step-size reduction
pars.rho     = 9.6*n; %different potential parameter (arbitrarily chosen)
pars.tol     = 1e-5;  %higher precision
[x,f,R]      = lpapgpot(A,b,c,pars);

fprintf('Results (high precision): \n');
fprintf('  %-25s%7.3e \n',' APGPOT Optimal value:',f);
fprintf('  %-25s%7.3e \n','LINPROG Optimal value:',f1);
