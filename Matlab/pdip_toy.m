clear all
format compact
%Sep 22 2012
%Primal dual interior point testbench. 
%This code is a simple implementation of the Mehorta 
%predictor corrector interior point method.

%It saves the generated linear systems and all 
%iteration data.

%It solves LPs in standard form
% min c'x s.t. Ax=b and 0 <= x

%{
%Uncomment this section to generate a random problem
%Problem definition
prob_name = 'Random';
seed = 0;
randn('seed',seed);

m = 100;
n = 1000;
A = randn(m,n);

x0    = rand(n,1);
b     = A*x0; %Find a rhs which makes the problem feasible

y0    = randn(m,1);
z0    = rand(n,1);
c     = A'*y0 + z0; %Generate a dual feasible
%}

%Uncomment this section to load problem from lpnetlib
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

x0 = 4000*ones(n,1);
y0 = zeros(m,1);
z0 = 4000*ones(n,1);

%Tolerances 
tol      = 1e-10;
max_iter = 50;

%Variables of the search direction
x_d = zeros(n,1);
y_d = zeros(m,1);
z_d = zeros(n,1);

%Should we use MINRES
iterative = false;
%minres options
shift = 0;
show  = false;
check = true;
itnlim= 200*m;
rtol  = 1e-16;


k  = 1; %iteration counter
go = true; %Iteration Flag

%Calculate the initial residuals
p_r  = b - A*x0;
d_r  = c - A'*y0 - z0;
c_r  = -x0.*z0;

%Calculate the initial centrality measure
mu   = x0'*z0;
mu   = mu/n;

n_p_r = norm(p_r);
n_d_r = norm(d_r);
n_c_r = norm(c_r);

%Log the step
log_x   = {x0};
log_y   = {y0};
log_z   = {z0};
log_mu  = {mu};
log_p_r = {p_r};
log_d_r = {d_r};
log_c_r = {c_r};

log_muaff = {};
log_sigma = {};

fprintf('%s, %s , %s,  %s,  %s,   %s,   %s, %s, %s,%s,%s,%s \n','Itn','itn_ass','itn_corr','n_p_r','n_d_r','mu','step x','step z','sigma','d_x','d_y','d_z');

%Iteration

%First iterates
x = x0;
y = y0;
z = z0;


while go
    %Form the 3 by 3 system
    z_h = sqrt(z);
    
    %Define the handle for the 3 by 3 system
    %K = @(d)([A'*d(n+1:n+m)+z_h.*d(n+m+1:2*n+m);A*d(1:n);z_h.*d(1:n)+x.*d(n+m+1:2*n+m)]);
     K = [[sparse(n,n),A',diag(z_h)];[A,sparse(m,n+m)];[diag(z_h),sparse(n,m),diag(x)]];


    %Define the RHS for the affine scaling problem
    RHS_aff = [d_r;p_r;-x.*z_h];
    if(iterative)
    %Solve the predictor system
        interrupt_save = @(varargin)({varargin{nargin}{:},varargin{nargin-1}});
        interrupt_minres = @(x,rnorm_est,interrupt_results)(interrupt_save(norm(RHS_aff-K*x),interrupt_results));

        [ x_sol, istop, itn_aff, rnorm, Arnorm, Anorm, Acond, ynorm, interrupt_results] = ...
                   minres( K, RHS_aff, [], shift, show, check, itnlim, rtol, interrupt_minres);
%         fprintf('Residual Norm Predictor %f, norm rhs %f, norm sol %f \n',log10(norm(K*x_sol-RHS_aff)),norm(RHS_aff),norm(x_sol));
    else
        itn_aff = 1;
        x_sol = K\RHS_aff; 
%        fprintf('Residual Norm Predictor %f, norm rhs %f, norm sol %f \n',log10(norm(K*x_sol-RHS_aff)),norm(RHS_aff),norm(x_sol));
    end

        clear RHS_aff;

        %Remove the symmetric scaling to x_sol
        x_sol(n+m+1:2*n+m) = z_h.*x_sol(n+m+1:2*n+m);
        
        x_aff = x_sol(1:n);
        y_aff = x_sol(n+1:n+m);
        z_aff = x_sol(n+m+1:2*n+m);

        %Find the largest afine step to the boundary
        max_alph_aff_x = -x./x_aff;
        max_alph_aff_x = min(max_alph_aff_x(find(max_alph_aff_x>0)));
        max_alph_aff_z = -z./z_aff;
        max_alph_aff_z = min(max_alph_aff_z(find(max_alph_aff_z>0)));
        
        %Find the predicted value of the duality measure
        mu_aff         = (x+max_alph_aff_x*x_aff)'*(z+max_alph_aff_z*z_aff); 
        mu_aff         = mu_aff/n;

        %Calculate the centering parameter 
        sigma = (mu_aff/mu).^3;
    
        %Solve the corrector system

        %Define the RHS for the corrector 
        z_m_h = 1./z_h;
        RHS_corr = [zeros(m+n,1);z_m_h.*(sigma*mu-x_aff.*z_aff)];
    if(iterative)
        %Solve the corrector system
        interrupt_save = @(varargin)({varargin{nargin}{:},varargin{nargin-1}});
        interrupt_minres = @(x,rnorm_est,interrupt_results)(interrupt_save(norm(RHS_corr-K*x),interrupt_results));

        [ x_sol, istop, itn_corr, rnorm, Arnorm, Anorm, Acond, ynorm, interrupt_results] = ...
                   minres( K, RHS_corr, [], shift, show, check, itnlim, rtol, interrupt_minres);
 %      fprintf('Residual Norm Corrector %f, norm solution %f, norm rhs %f \n',log10(norm(K*x_sol-RHS_corr)),log10(norm(x_sol)),log10(norm(RHS_corr)));
    else
        itn_corr = 1;
        x_sol =K\RHS_corr;         
 %       fprintf('Residual Norm Corrector %f, norm solution %f, norm rhs %f \n',log10(norm(K*x_sol-RHS_corr)),log10(norm(x_sol)),log10(norm(RHS_corr)));
    end     
        clear RHS_corr;

        %Remove the symmetric scaling to x_sol
        x_sol(n+m+1:2*n+m) = z_h.*x_sol(n+m+1:2*n+m);
        
        x_corr = x_sol(1:n);
        y_corr = x_sol(n+1:n+m);
        z_corr = x_sol(n+m+1:2*n+m);

%        fprintf('Axcorr: %f \nAtycorr + zcorr - dr: %f \nXzcorr + Zxcorr +Xzcorr %f\n',norm(A*x_corr),norm(A'*y_corr+z_corr),norm(z.*x_corr + x.*z_corr + x_aff.*z_aff - sigma*mu));
    %Sanity check, if the linear system is solved accurately, the x_corr and y_corr
    %directions must be close to the null space of their equations
    
    %Calculate the directions
    x_d = x_aff+x_corr;
    y_d = y_aff+y_corr;
    z_d = z_aff+z_corr;

    %Find the largest afine step to the boundary
    max_alph_x = -x./x_d;
    max_alph_x = min(max_alph_x(find(max_alph_x>0)));
    max_alph_z = -z./z_d;
    max_alph_z = min(max_alph_z(find(max_alph_z>0)));
    
    %Calculate the step length
    alph_x = min(0.99*max_alph_x,1);
    alph_z = min(0.99*max_alph_z,1); 
    
    if(length(alph_z)==0)
        fprintf('______ problem seems to be unbounded alph_z == 0--------');
        return
    elseif(length(alph_x) == 0)
        fprintf('______problem seems unbounded alph_x == 0-----'); 
        return
    end
    %update the variables
    x  = x+alph_x*x_d;
    y  = y+alph_z*y_d;
    z  = z+alph_z*z_d;
  
    %Calculate the new residuals
    p_r  = b - A*x;
    d_r  = c - A'*y - z;
    c_r  = -x.*z;

    %Calculate the new centrality measure
    mu   = x'*z;
    mu   = mu/n;
  
    %Stopping criteria 
    n_p_r = norm(p_r);
    n_d_r = norm(d_r);
    n_c_r = norm(c_r);
   
    %Log the step
    log_x   = {log_x{:},x};
    log_y   = {log_y{:},y};
    log_z   = {log_z{:},z};
    log_mu  = {log_mu{:},mu};
    log_p_r = {log_p_r{:},p_r};
    log_d_r = {log_d_r{:},d_r};
    log_c_r = {log_c_r{:},c_r};
    log_muaff = {log_muaff{:},mu_aff};
    log_sigma = {log_sigma{:},sigma};

    %Residual sum
    res_sum = n_p_r + n_d_r + n_c_r;

    %Print some stuff
    fprintf('%3i, %7i ,%7i , %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f \n',...
        k,itn_aff,itn_corr,log10(n_p_r),log10(n_d_r),log10(mu),alph_x,alph_z,sigma,norm(x_d),norm(y_d),norm(z_d));

    if res_sum < tol || k == max_iter
        go = false;
    else
        k = k+1;
    end
    
end
%remove the backslashes from the name
prob_name(find(prob_name=='/'))='_';
prob_name = ['Runs/',prob_name];
save(prob_name,'k','log_d_r','log_p_r','log_mu','log_x','log_y','log_z','log_muaff','log_sigma','A','b','c','prob_name');
