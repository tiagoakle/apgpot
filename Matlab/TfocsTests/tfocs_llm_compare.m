%Executes Lan, Lu, Monteiros aproach using TFOCS

%Jan 29 2013

%This file assumes that the problem is defined
%in the variables A,b,c

%Define the objective min||Es+q||
m = size(A,1);
n = size(A,2);

E = [A      ,zeros(m,m),zeros(m,n) ;...
    zeros(n),A'        ,eye(n)     ;...
    c'      ,-b'       ,zeros(1,n)];
q = [-b;-c;0];

%Lan lu monteiro criterion
pweight = max(1,norm(b));
dweight = max(1,norm(c));
stoptol = 0.01;
err_obj             = @(f,x)(c'*x(1:n)); %function to store the objective value of the lp
err_primal_rel_inf  = @(f,x)(norm(A*x(1:n)-b)/pweight);
err_dual_rel_inf    = @(f,x)(norm(A'*x(n+1:n+m) + x(n+m+1:2*n+m) - c)/dweight); 
err_rel_gap         = @(f,x)((c'*x(1:n)-b'*x(n+1:n+m))/(max(1,0.5*abs(c'*x(1:n))+0.5*abs(b'*x(n+1:n+m)))));
stop_crit           =@(f,x)((err_primal_rel_inf(f,x)<stoptol)&&(err_rel_gap(f,x)<stoptol)&&(err_dual_rel_inf(f,x)<stoptol));


opts.saveHist = true;
%opts.alg      = 'LLM';
opts.errFcn   = {err_obj,err_primal_rel_inf,err_dual_rel_inf,err_rel_gap};
opts.stopFcn  = stop_crit;

[sol,res] = tfocs(smooth_quad,{E,q},proj_restrictedRplus(n,m),ones(2*n+m,1),opts);


