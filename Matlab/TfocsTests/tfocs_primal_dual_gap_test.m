%Generates a problem in standard form and tries to solve it 
%minimizing a constrained quadratic in TFOCS.
%Jan 25 2013

%If there is a problem in the workspace then dont generate a new one
if(~exist('A','var')||~exist('b','var')||~exist('c','var'))
    %Problem definition
    prob_name = 'Random';
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
end
%Problem parameters
m = size(A,1);
n = size(A,2);

E = [A      ,zeros(m,m),zeros(m,n) ;...
    zeros(n),A'        ,eye(n)     ;...
    c'      ,-b'       ,zeros(1,n)];

opts.saveHist = true;
opts.maxIts   = 1000;
opts.alg      = 'AT'
opts.errFcn   = @(f,x)(c'*x(1:n)); %function to store the objective value of the lp
[sol,res,out_pars] = tfocs(smooth_quad,{E},proj_restrictedRplus(n,m),ones(2*n+m,1),opts);

plot(res.err);
xsol = sol(1:n);
ysol = sol(n+1:n+m);
zsol = sol(n+m+1:2*n+m);

%Re solve using the lp solver
obj    = smooth_linear(c);

x0     = [];
z0     = [];
mu     = 1.e-5;
opts   = struct;
opts.saveHist = true;
opts.alg      = 'AT'

[sol,res] = solver_sLP( c, A,-b, mu);

