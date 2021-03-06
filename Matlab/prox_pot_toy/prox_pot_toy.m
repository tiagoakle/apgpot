%clear all
%clc
format compact

%December 5 2012
%Primal dual proximal gradient mapping potential reduction method
%This code implements the experimental code from Zhuo, Skaja, Ye, Akle. 

%Jan 16 2013
%The code now has options to choose between the two objective functions f1,f2 and
%to choose between projecting a la potential reduction or interpreting h as an indicator function

%It saves the generated gradients and all
%iteration data.

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
if~exist('pars','var')
  pars = struct;
end


% set defaults:
defaults = {...
    'accel',1;... %Use nesterov acceleration
    'epsi',1e-10;... %Tolerance to stop
    'max_bkt', 100;... %Backtracking iterations
    'beta_bkt',0.5;... %Backtracking constant
    't0',1;...         %Initial step size
    't0m',1/0.9;...    %Extension to the step size
    'func','gap';...   %Function to use to solve lp [gap,cent]
    'proj','pot';...   %Projection type [pot,ind]
    'restart',0;...    %Restart momentum every so many iterations 0 is no restart
    'max_iter',1e6;... %Iteration limit
    'itn_print',1000;  %Print every so many iterations
    'rho',-1;          %Rho -1 means that rho will be chosen as n+sqrt(n)
 };

for k = 1:size(defaults,1)
    if ~isfield(pars,defaults(k,1))
        pars.(defaults{k,1}) = defaults{k,2};
    end
end

%Check for some invalid parameters
if(strcmp(pars.func,'phi')&&strcmp(pars.proj,'ind'))
    error('ProxPotToy:Invalid choice of projection for objective phi');
end

%Define the variables for the log
results.log_mu        = zeros(pars.max_iter,1);
results.log_p_r       = zeros(pars.max_iter,1);
results.log_d_r       = zeros(pars.max_iter,1);
results.log_g_r       = zeros(pars.max_iter,1);
results.log_merit     = zeros(pars.max_iter,1);
results.log_f         = zeros(pars.max_iter,1);
results.log_step      = zeros(pars.max_iter,1);
results.log_infeas    = zeros(pars.max_iter,1);
results.log_rel_infeas= zeros(pars.max_iter,1);

results.exit_flag = 'max_iter';

if pars.rho == -1
    rho = (n+sqrt(n)); %Feasibility weight in the potential function
else
    rho = pars.rho;
end

itn  = 1; %iteration counter
k    = 0; %Backtraking interation couter
rescount       = 1; %Counter since last momentum restart

%Calculate the initial residuals
p  = A*x0-b;
d  = A'*y0 + z0 - c;
g  = c'*x0-b'*y0;
o  = x0'*c;
ce = (x0'*z0);
mu = ce/n;  


%Measure the infeasibilities
nP0 = norm(p);
nD0 = norm(d);
nG0 = abs(g);

nP  = nP0;
nD  = nD0;
nG  = nG0;

res0   = norm([p;d;g]);

%Calculate the initial value of f
if strcmp(pars.func,'gap')
    f       = nP^2+nD^2+nG^2;   
elseif strcmp(pars.func,'cent')
    f       = nP^2+nD^2+ce;    
end


%First iterates
x = x0;
y = y0;
z = z0;
xp= x;
yp= y;
zp= z;

f0 = f;
      
%First step size
t  = pars.t0;

fprintf('Proximal potential reduction toy code\n')
fprintf('LP objective function %s \n',pars.func);
fprintf('Projection into positive orthant %s \n',pars.proj)
fprintf('Acceleration %i \n',pars.accel);
fprintf('Rho: %d\n',pars.rho);

fprintf('================================ PROXPOT TOY =========================================================================\n');
fprintf('  ITE LSTP F(X)/F0  f(x)/f(0) PHI     Mu       NPR      NDR      NGR      ALPH     MinX     MinZ     Infeas   R_Infeas\n');
fprintf('----------------------------------------------------------------------------------------------------------------------\n');

%-----------------------------------------------------
% Start the iteration
%----------------------------------------------------

for itn = 1:pars.max_iter

   %Exit condition:
    if nP < pars.epsi && nD < pars.epsi && mu < pars.epsi
     fprintf('Success \n');
     results.exit_flag = 'Success';
     break;
    end
    
    %----------------------------- 
    %Evaluate the progress and log
    %----------------------------- 
    
    relf  = f/f0; %Relative reduction in objective f
    if strcmp(pars.func,'gap') 
        merit = rho*log(nP^2 + nD^2 + nG^2) - sum(log(x)) - sum(log(z)); %Value of phi
    elseif strcmp(pars.func,'cent') 
        merit = rho*log(nP^2 + nD^2 + ce) - sum(log(x)) - sum(log(z)); %Value of phi
    end


    res   = norm([p;d;g]); %Norm of the residual
    relres= res/res0;      %Relative residual 

    %calculate the objetive value
    obj   = x'*c;
    %Calculate the centrality 
    ce    = z'*x;
    mu    = ce/n;
    %Calculate the relative infeasibility
    rel_infeas = max([nP/nP0,nD/nD0,nG/nG0]); %Relative infeasibility/ to the original
    infeas = max([nP,nD,nG]);

    %Log the step
    results.log_mu(itn)        = mu;
    results.log_p_r(itn)       = nP;
    results.log_d_r(itn)       = nD;
    results.log_g_r(itn)       = nG;
    results.log_merit(itn)     = merit;
    results.log_f(itn)         = f;
    results.log_step(itn)      = t;
    results.log_infeas(itn)    = infeas;
    results.log_rel_infeas(itn)= rel_infeas;
    results.log_obj(itn)       = o;

   if mod(itn,pars.itn_print) == 0 || itn ==1
            fprintf('%5.1i %3.1i %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e\n',...
                     itn,k,relres,relf,merit,mu,nP/nP0,nD/nD0,nG/nG0,t,min(x),min(z),infeas,rel_infeas); 
   end

   
    %---------------------- 
    %Form the gradient of f
    %----------------------
    if strcmp(pars.func,'gap')
        grad = 2*[A'*p+g*c;A*d-g*b;d]; 
    elseif strcmp(pars.func,'cent')
        grad = 2*[A'*p;A*d;d] + [z;zeros(m,1);x];
    end

    
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
        if(pars.proj == 'pot') %Potential reduction projection
            eta = 4*t*f/rho;
            x   = 0.5*(x+sqrt(x.^2+eta)); 
            z   = 0.5*(z+sqrt(z.^2+eta));
        elseif(pars.proj == 'ind') %Positive orthant projection
            x   = max(x,0);
            z   = max(z,0); 
        end
                
        %Calculate the residuals
        p = A*x-b;
        d = A'*y+z-c;
        g = (c'*x-b'*y);
        o = x'*c;
        ce= x'*z;
        mu= ce/n;

        %Calculate the infeasibility 
        nP = norm(p);
        nD = norm(d);
        nG = abs(g);
        
        %Calculate the proximal gradient
        Gt  = 1/t*[xs-x;ys-y;zs-z];
        Gtg = Gt'*grad; %local linear decrease
        nGt = norm(Gt)^2;


       if strcmp(pars.func,'gap')
            fn  = nP^2+nD^2+g^2;  
       elseif strcmp(pars.func,'cent')
            fn  = nP^2+nD^2+ce;  
       end

        
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
        g = (c'*x-b'*y);
        o = x'*c;
        ce= x'*z;
        mu= ce/n;
        %Calculate the infeasibility 
        nP = norm(p);
        nD = norm(d);
        nG = abs(g);
           
        if strcmp(pars.func,'gap')
            f  = nP^2+nD^2+g^2;  
        elseif strcmp(pars.func,'cent')
            f  = nP^2+nD^2+ce;  
        end

    end

 
    %-------------------------------------
    %Fixed restarting 
    %-------------------------------------
    if  pars.restart > 0 && pars.restart == rescount
      xp = x;
      yp = y;
      zp = z;
      'Restart'
      rescount = 0;
    else
      rescount = rescount +1;
    end

   %Prepare the next loop
   if(isfield(pars,'t0m')) %Extend the initial step size a little
        t = t*pars.t0m;
   end
          
end

%---------------------------------------------------------
%End of the iteration
%---------------------------------------------------------
results.itn = itn;
