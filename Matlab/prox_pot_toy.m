clear all
format compact

%December 5 2012
%Primal dual proximal gradient mapping potential reduction method
%This code implements the experimental code from Zhuo, Skaja, Ye, Akle. 

%It saves the generated gradients and all
%iteration data.

%It solves LPs in standard form
% min c'x s.t. Ax=b and 0 <= x

%Problem definition
prob_name = 'Random';
seed = 0;
randn('seed',seed);

m = 200;
n = 1000;
A = randn(m,n);

x0    = rand(n,1);
b     = A*x0; %Find a rhs which makes the problem feasible

y0    = randn(m,1);
z0    = rand(n,1);
c     = A'*y0 + z0; %Generate a dual feasible

%Problem parameters
m = size(A,1);
n = size(A,2);

%{
%Load the indices of the standard form 
%problems in lpnetlib
load 'standard_form_indices';
problem_index = 2;
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
%}

x0 = ones(n,1);
y0 = zeros(m,1);
z0 = ones(n,1);

%Tolerances 
epsi      = 1e-10;
max_iter  = 3000;

%Variables of the search direction
x_d = zeros(n,1);
y_d = zeros(m,1);
z_d = zeros(n,1);

rho = (n+sqrt(n)); %Feasibility weight in the potential function


itn  = 1; %iteration counter
go = true; %Iteration Flag

%Calculate the initial centrality measure
mu   = x0'*z0;
mu   = mu/n;

%Calculate the initial residuals
p  = A*x0-b;
d  = A'*y0 + z0 - c;
g  = c'*x0-b'*y0;
mu = (x0'*z0)/(n);  
%Calculate the infeasibility 
f = norm(p)^2+norm(d)^2+g^2;    
merit = rho*log(norm(p)^2 + norm(d)^2 + g^2) - sum(log(x0)) - sum(log(z0));
%----------------------
%Tuning variables
%----------------------
stepeps = 1e-5;
lambda  = 1;

%Log the step
log_x   = {x0};
log_y   = {y0};
log_z   = {z0};
log_mu  = {mu};
log_p_r = {p};
log_d_r = {d};
log_g_r = {g};
log_merit = {merit};
log_f   = {f};
log_step= {};
fprintf('%s, %s , %s,  %s,  %s,   %s,   %s, %s, %s,%s,%s,%s \n','Itn','itn_ass','itn_corr','n_p_r','n_d_r','mu','step x','step z','sigma','d_x','d_y','d_z');

%Iteration

%First iterates
x = x0;
y = y0;
z = z0;


while go

    %Form the hsd system 
    % K = [[sparse(n,n),c,-A',-speye(n),sparse(n,1)];[-c',0,b',sparse(1,n),-1];[A,b,sparse(m,n+m+1)];[diag(sparse([z;k])),sparse(n+1,m),diag(sparse([x;t]))]];

    %Define the RHS for the affine scaling problem
    %RHS = [eta*d;eta*g;-eta*p;sigma*mu-x.*z;sigma*mu-t*k];
    %x_sol = K\RHS; 
    %clear RHS;
   
    %----------------- 
    %Form the gradient
    %-----------------
    grad = rho/f*[A'*p+g*c;A*d-g*b;d];

    
    %Calculate the local decrease
    dec  = grad'*[x;y;z];

    %Save the present step    
    xs = x;
    ys = y;
    zs = z;

    %Calculate the step size
    %XXX
    alph = stepeps/(2*rho*lambda);

    %Very aggressive strategy -------------------------

    %Find the largest afine step to the boundary
   % max_alph_x = -[x]./[grad(1:n)];
   % max_alph_x = min(max_alph_x(find(max_alph_x>0)));
   % max_alph_z = -[z]./[grad(n+m+1:2*n+m)];
   % max_alph_z = min(max_alph_z(find(max_alph_z>0)));
   % max_alph   = min([max_alph_z,max_alph_x,1]);
 
   % %Calculate the step length
   % alph_x = min(0.90*max_alph_x,1);
   % alph_z = min(0.90*max_alph_z,1); 
   % %Choose the minimum last feasible step
   % alph   = min(alph_x,alph_z);
   % alph   = min(alph,1);
  

    %-------------------------------
    %Take the gradient step
    %-------------------------------
    x  = x - alph*grad(1:n);
    y  = y - alph*grad(n+1:n+m);
    z  = z - alph*grad(n+m+1:2*n+m); 

    %--------------------------------
    %Proximal projection step
    %--------------------------------
    x  = 0.5*(x+sqrt(x.^2+4*alph));
    y  = 0.5*(y+sqrt(y.^2+4*alph));
    z  = 0.5*(z+sqrt(z.^2+4*alph));

    %Calculate the residuals
    p = A*x-b;
    d = A'*y+z-c;
    g = (c'*x-b'*y);
    %Calculate the infeasibility 
    f = norm(p)^2+norm(d)^2+g^2;  

    pastmerit = merit;
    merit = rho*log(norm(p)^2 + norm(d)^2 + g^2) - sum(log(x)) - sum(log(z));
    %Calculate the centrality 
    mu = x'*z/n;

    %Log the step
    log_x     = {log_x{:},x};
    log_y     = {log_y{:},y};
    log_z     = {log_z{:},z};
    log_mu    = {log_mu{:},mu};
    log_p_r   = {log_p_r{:},p};
    log_d_r   = {log_d_r{:},d};
    log_g_r   = {log_g_r{:},g};
    log_merit = {log_merit{:},merit};
    log_f     = {log_f{:},f};
    log_step  = {log_step{:},alph};
     
   fprintf(' merit %d ,dec %d, ini-final %d, step %d, f %d , x (%d,%d), z(%d,%d) \n',merit,dec,merit-pastmerit,alph,f,min(x),max(x),min(z),max(z));

  %Update mu
  mu = (x'*z)/n;
  if itn > max_iter 
      go = false;
  else
      itn = itn +1;
  end
  if f < epsi
    go = false;
    fprintf('Success \n');
  end

end
%remove the backslashes from the name
prob_name(find(prob_name=='/'))='_';
%save(prob_name,'k','log_d_r','log_p_r','log_mu','log_x','log_y','log_z','log_muaff','log_sigma');

