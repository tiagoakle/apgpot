%clear all
%clc
format compact

%December 5 2012
%Primal dual proximal gradient mapping potential reduction method
%This code implements the experimental code from Zhuo, Skaja, Ye, Akle. 

%It saves the generated gradients and all
%iteration data.

%It solves LPs in standard form
% min c'x s.t. Ax=b and 0 <= x

%If there is a problem in the workspace then dont generate a new one
if(isempty(A)||isempty(b)||isempty(c))
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
y0 = ones(m,1);
z0 = ones(n,1);

%---------------------------------------
%Define the parameters
%---------------------------------------
pars.accel     = 1;  %Use Nesterov acceleration 
%Tolerances 
pars.epsi      = 1e-5;
%backtracking
pars.max_bkt   = 100;
pars.beta_bkt  = 0.5;
pars.t0        = 1;  %Inital step size
pars.t0m       = 1/0.9; %Extension of the step size

%printing
Exit_flag = 'max_iter';
rho = (n+sqrt(n)); %Feasibility weight in the potential function

itn  = 1; %iteration counter
k    = 0; %Backtraking interation couter
%Calculate the initial centrality measure
mu   = x0'*z0;
mu   = mu/n;

%Calculate the initial residuals
p  = A*x0-b;
d  = A'*y0 + z0 - c;
g  = sqrt(x.*z);
o  = x0'*c;

mu = (x0'*z0)/(n);  
%Measure the infeasibilities
nP0 = norm(p);
nD0 = norm(d);
nG0 = norm(g);

nP  = nP0;
nD  = nD0;
nG  = nG0;

res0   = norm([p;d;g]);
%Calculate the infeasibility 
f       = nP^2+nD^2+nG^2;    

%Log the step
log_mu    = zeros(pars.max_iter,1);
log_p_r   = zeros(pars.max_iter,1);
log_d_r   = zeros(pars.max_iter,1);
log_g_r   = zeros(pars.max_iter,1);
log_merit = zeros(pars.max_iter,1);
log_f     = zeros(pars.max_iter,1);
log_step  = zeros(pars.max_iter,1);
log_infeas= zeros(pars.max_iter,1);
log_rel_infeas= zeros(pars.max_iter,1);
log_obj       = zeros(pars.max_iter,1);

%Iteration

%First iterates
x = x0;
y = y0;
z = z0;
xp= x;
yp= y;
zp= z;

f0 = f;
      
%XXX: the choice is arbitrary
%First step size
t  = 1;

fprintf('================================ XZ PROXPOT TOY ======================================================================\n');
fprintf('  ITE LSTP F(X)/F0  f(x)/f(0) PHI     Mu       NPR      NDR      NGR      ALPH     MinX     MinZ     Infeas   R_Infeas\n');
fprintf('----------------------------------------------------------------------------------------------------------------------\n');

k       = 1; %Counter since momentum restart
for itn = 1:pars.max_iter

   %Exit condition:
    if nP < pars.epsi && nD < pars.epsi && mu < pars.epsi
     fprintf('Success \n');
     exit_flag = 'Success';
     break;
    end
    
    %----------------------------- 
    %Evaluate the progress and log
    %----------------------------- 
    
    relf = f/f0; %Relative reduction in objective f
    merit = rho*log(nP^2 + nD^2 + nG^2) - sum(log(x)) - sum(log(z)); %Value of phi
    res   = norm([p;d;g]); %Norm of the residual
    relres= res/res0;      %Relative residual
    %calculate the objetive value
    obj= x'*c;
    %Calculate the centrality 
    mu = x'*z/n;
    rel_infeas = max([nP/nP0,nD/nD0,nG/nG0]); %Relative infeasibility/ to the original
    infeas = max([nP,nD,nG]);

    %Log the step
    log_mu(itn)        = mu;
    log_p_r(itn)       = nP;
    log_d_r(itn)       = nD;
    log_g_r(itn)       = nG;
    log_merit(itn)     = merit;
    log_f(itn)         = f;
    log_step(itn)      = t;
    log_infeas(itn)    = infeas;
    log_rel_infeas(itn)= rel_infeas;
    log_obj(itn)       = o;

   if mod(itn,pars.itn_print) == 0 || itn ==1
            fprintf('%5.1i %3.1i %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e\n',...
                     itn,k,relres,relf,merit,mu,nP/nP0,nD/nD0,nG/nG0,t,min(x),min(z),infeas,rel_infeas); 
   end

   
    %---------------------- 
    %Form the gradient of f
    %----------------------
    grad = 2*[A'*p;A*d;d] + [z;zeros(m,1);x];
    
    %Save the present step    
    xs = x;
    ys = y;
    zs = z;
    
    %-----------------------------------
    %Backtracking linesearch
    %----------------------------------

    for k = 1:pars.max_bkt
        %-------------------------------
        %Take the gradient step
        %-------------------------------
        x  = xs - t*grad(1:n);
        y  = ys - t*grad(n+1:n+m);
        z  = zs - t*grad(n+m+1:2*n+m); 

        %--------------------------------
        %Proximal projection step on the 
        % positively constrained variables
        %-------------------------------
        %XXX This is to win intuition for the proof
        %eta = 4*t*1.e-10/rho; %OOPS this works too!
        
        eta = 4*t*f/rho;
        x   = 0.5*(x+sqrt(x.^2+eta)); 
        z   = 0.5*(z+sqrt(z.^2+eta));
        
        %Calculate the residuals
        p = A*x-b;
        d = A'*y+z-c;
        g = sqrt(x.*z);
        o = x'*c;

        %Calculate the infeasibility 
        nP = norm(p);
        nD = norm(d);
        nG = norm(g);

        %Calculate the proximal gradient
        Gt  = 1/t*[xs-x;ys-y;zs-z];
        Gtg = Gt'*grad; %local linear decrease
        nGt = norm(Gt)^2;
        fn  = nP^2+nD^2+nG^2;  
        
        %Evaluate the backtracking inequality 
        if fn < f - t*Gtg + t*0.5*nGt
             break;
        else %Backtracking step
            t = t*pars.beta_bkt;
        end

    end % end of linesearch

    f = fn;

    %---------------------------------
    %Nesterov/ P.Tseng Acceleration
    %--------------------------------- 
    if pars.accel > 0  %Acceleration step
        theta   = (itn-1)/(itn+2); %Extrapolation coeficient
        ex_dirx = [x-xp];
        ex_diry = [y-yp];
        ex_dirz = [z-zp];
       
        %Save x
        xp      = x;
        zp      = z;
        yp      = y;
        resp    = res;
        %Extrapolate
        x     = x + theta*ex_dirx;
        y     = y + theta*ex_diry;
        z     = z + theta*ex_dirz; 
        %Update the residuals
        %XXX:The expensive way
        p = A*x-b;
        d = A'*y+z-c;
        g = sqrt(x.*z);
        o = x'*c;
        
        %Calculate the infeasibility 
        nP = norm(p);
        nD = norm(d);
        nG = norm(g);

        f  = nP^2+nD^2+nG^2;  
        
    end
    
    %-------------------------------------
    %Fixed restarting 
    %-------------------------------------
    if isfield(pars,'restart') && mod(itn,pars.restart) == 0
        xp = x;
        yp = y;
        zp = z;
        k = 0;
        'Restart'
    end

   %Prepare the next loop
   if(isfield(pars,'t0m')) %Extend the initial step size a little
        t = t*pars.t0m;
   end
          
k = k+1;
end
%remove the backslashes from the name
%prob_name(find(prob_name=='/'))='_';
%save(prob_name,'k','log_d_r','log_p_r','log_mu','log_x','log_y','log_z','log_muaff','log_sigma');

