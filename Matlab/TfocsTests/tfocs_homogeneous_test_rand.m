%Jan 22 2013
%Use tfocs to solve the homogeneous self dual model

%Jan 29 Modified to solve an lpnetlib example, while storing the objective value 
%   primal and dual infeasibility, cetnrality, gap, tau and kappa per iteration.
%It solves LPs in standard form
% min c'x s.t. Ax=b and 0 <= x

%Feb 01 Compare the behavior of HSD between a random an lpnetlib problem. The 
%sparsity of the random problem is taken from the lpnetlib one.

close all;
clear all
%Load an LPNETLIB problem
P = UFget(622);
A_lpnetlib = P.A;
c_lpnetlib = P.aux.c;
b_lpnetlib = P.b;

%%If there is a problem in the workspace then dont generate a new one
if(~exist('A','var')||~exist('b','var')||~exist('c','var'))
    %Problem definition
    prob_name = 'Random';
    seed = 0;
    randn('seed',seed);
    rand('seed',seed);
    m = 30;
    n = 100;
    density = 100;
    A = sprand(A_lpnetlib);

    x0    = rand(n,1);
    b     = A*x0; %Find a rhs which makes the problem feasible

    y0    = randn(m,1);
    z0    = rand(n,1) + 1;
    c     = A'*y0 + z0; %Generate a dual feasible
end

%Problem parameters
m = size(A,1);
n = size(A,2);
 
E = [A      ,-b,zeros(m,m),zeros(m,n),zeros(m,1);...
    zeros(n),-c,A'        ,eye(n)    ,zeros(n,1);...
    c'      ,0 ,-b'       ,zeros(1,n),1        ];

x_start = ones(2*(n+1)+m,1);

opts.saveHist = true;
opts.maxIts   = 10000;
opts.alg      = 'LLM'
opts.tol      = 0;%Go forever

err_obj       = @(f,x)(c'*x(1:n)); %function to store the objective value of the lp
err_primalIn  = @(f,x)(norm(A*x(1:n)-x(n+1)*b));
err_dualIn    = @(f,x)(norm(A'*x(n+2:n+m+1) + x(n+m+2:2*n+m+1) - x(n+1)*c)); 
err_central   = @(f,x)(x(1:n+1)'*x(n+m+2:2*n+m+2));
err_gap       = @(f,x)(c'*x(1:n)-b'*x(n+2:n+m+1));
err_tau       = @(f,x)(x(n+1));
err_kappa     = @(f,x)(x(2*(n+1)+m));
opts.errFcn   = {err_obj,err_primalIn,err_dualIn,err_central,err_gap,err_tau,err_kappa};
[sol_1,res] = tfocs(smooth_quad,{E},proj_restrictedRplus(n +1,m),x_start,opts);
sol_x = sol_1(1:n);
sol_t = sol_1(n+1);
sol_y = sol_1(n+2:n+m+1);
sol_z = sol_1(n+m+2:2*n+m+1);
sol_k = sol_1(2*n+m+2);

figure;
semilogy(abs(res.err(:,1)));
hold on
legend('Objective');
hold off

h = semilogy(abs(res.err(:,2)));
hold on
semilogy(abs(res.err(:,3)),'r');
legend('Primal','Dual');
saveas(h,'../../Writeup/Images/hsd_llm_primal_dual_infeas.png');
hold off

tauoverkappa = res.err(:,6)./res.err(:,7);
figure
h = semilogy(tauoverkappa);
hold on
legend('\tau \kappa');
hold off
saveas(h,'../../Writeup/Images/hsd_llm_tau_over_kappa.png');

figure;
h=semilogy(abs(res.err(:,4)));
hold on
legend('Centrality');
hold off
saveas(h,'../../Writeup/Images/hsd_llm_centrality.png');

figure;
semilogy(abs(res.err(:,5)));
hold on
legend('Gap');
hold off

%--------------REPEAT THE EXPERIMENT WITH THE LPNETLIB PROBLEM

A = A_lpnetlib;
b = b_lpnetlib;
c = c_lpnetlib;

E = [A      ,-b,zeros(m,m),zeros(m,n),zeros(m,1);...
    zeros(n),-c,A'        ,eye(n)    ,zeros(n,1);...
    c'      ,0 ,-b'       ,zeros(1,n),1        ];

x_start = ones(2*(n+1)+m,1);

opts.saveHist = true;
opts.maxIts   = 10000;
opts.alg      = 'LLM'
opts.tol      = 0;%Go forever

err_obj       = @(f,x)(c'*x(1:n)); %function to store the objective value of the lp
err_primalIn  = @(f,x)(norm(A*x(1:n)-x(n+1)*b));
err_dualIn    = @(f,x)(norm(A'*x(n+2:n+m+1) + x(n+m+2:2*n+m+1) - x(n+1)*c)); 
err_central   = @(f,x)(x(1:n+1)'*x(n+m+2:2*n+m+2));
err_gap       = @(f,x)(c'*x(1:n)-b'*x(n+2:n+m+1));
err_tau       = @(f,x)(x(n+1));
err_kappa     = @(f,x)(x(2*(n+1)+m));
opts.errFcn   = {err_obj,err_primalIn,err_dualIn,err_central,err_gap,err_tau,err_kappa};
[sol_1,res] = tfocs(smooth_quad,{E},proj_restrictedRplus(n +1,m),x_start,opts);
sol_x = sol_1(1:n);
sol_t = sol_1(n+1);
sol_y = sol_1(n+2:n+m+1);
sol_z = sol_1(n+m+2:2*n+m+1);
sol_k = sol_1(2*n+m+2);

figure;
semilogy(abs(res.err(:,1)));
hold on
legend('Objective');
hold off

h = semilogy(abs(res.err(:,2)));
hold on
semilogy(abs(res.err(:,3)),'r');
legend('Primal','Dual');
saveas(h,'../../Writeup/Images/hsd_llm_primal_dual_infeas_lpentlib.png');
hold off

tauoverkappa = res.err(:,6)./res.err(:,7);
figure
h = semilogy(tauoverkappa);
hold on
legend('\tau \kappa');
hold off
saveas(h,'../../Writeup/Images/hsd_llm_tau_over_kappa_lpentlib.png');

figure;
h=semilogy(abs(res.err(:,4)));
hold on
legend('Centrality');
hold off
saveas(h,'../../Writeup/Images/hsd_llm_centrality_lpentlib.png');

figure;
semilogy(abs(res.err(:,5)));
hold on
legend('Gap');
hold off


