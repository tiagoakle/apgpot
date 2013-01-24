%Jan 22 2013
%Use tfocs to solve the homogeneous self dual model

%It solves LPs in standard form
% min c'x s.t. Ax=b and 0 <= x

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

E = [A      ,-b,zeros(m,m),zeros(m,n),zeros(m,1);...
    zeros(n),-c,A'        ,eye(n)    ,zeros(n,1);...
    c'      ,0 ,-b'       ,zeros(1,n),1        ];

opts.saveHist = true;
opts.maxIts   = 1000;
opts.alg      = 'AT'
[sol,res] = tfocs(smooth_quad,{E},proj_restrictedRplus(n+1,m),ones(2*(n+1)+m,1),opts);
semilogy(res.f);
xsol = sol(1:n);
tsol = sol(n+1);
ysol = sol(n+2:n+m+1);
zsol = sol(n+m+2:2*n+m+1);
ksol = sol(2*n+2+m);
err = E*sol;
nerr = norm(err);
