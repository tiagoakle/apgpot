%clear all
%clc
format compact

%Feb 1 2013
%This code does accelerated gradient descent without projection

%It saves the generated gradients and all
%iteration data.

%It solves LPs in standard form
% min c'x s.t. Ax=b and 0 <= x

%Define problem parameters
m=524;
n=1028;
density = 2;
%m = 5;
%n = 15;
%density = 100;

%Generate the problem instance
A = sprand(m,n,density/100);
A = randn(m,n);
x0    = rand(n,1);
b     = A*x0; %Find a rhs which makes the problem feasible 
y0    = randn(m,1);
z0    = rand(n,1);
c     = A'*y0 + z0; %Generate a feasible dual

%
%
%%If there is a problem in the workspace then dont generate a new one
%if(~exist('A','var')||~exist('b','var')||~exist('c','var'))
%    %Problem definition
%    prob_name = 'Random';
%    seed = 0;
%    randn('seed',seed);
%    rand('seed',seed);
%    m = 50;
%    n = 150;
%    A = randn(m,n);
%
%    x0    = rand(n,1);
%    b     = A*x0; %Find a rhs which makes the problem feasible
%
%    y0    = randn(m,1);
%    z0    = rand(n,1);
%    c     = A'*y0 + z0; %Generate a dual feasible
%end
%%Problem parameters
%m = size(A,1);
%n = size(A,2);

%
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
    't0m',1/0.6;...    %Extension to the step size
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
f   = nP^2+nD^2+nG^2;
ps  = 2*rho*log(f) - sum(log(x0)) - sum(log(z0));

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

if(pars.itn_print>0)
    fprintf('Proximal potential reduction by experimental method 2\n')
    fprintf('Acceleration %i \n',pars.accel);
    fprintf('Rho: %d\n',pars.rho);
    
    
    fprintf('================================ PROXPOT TOY Experimental ============================================================\n');
    fprintf('  ITE LSTP F(X)/F0  f(x)/f(0) PHI     Mu       NPR      NDR      NGR      ALPH     MinX     MinZ     Infeas   R_Infeas\n');
    fprintf('----------------------------------------------------------------------------------------------------------------------\n');
end

%-----------------------------------------------------
% Start the iteration
%----------------------------------------------------
stop_in_error = false;
for itn = 1:pars.max_iter

   %Exit condition:
    if nP < pars.epsi && nD < pars.epsi && mu < pars.epsi
     if(pars.itn_print > 0)
        fprintf('Success \n');
     end
     results.exit_flag = 'Success';
     break;
    end
    if(isfield(pars,'stopFnc'))
        if(pars.stopFnc(x,y,z))
          results.exit_flag = 'User specified exit criteria';
          break;
        end

    end

    
    %----------------------------- 
    %Evaluate the progress and log
    %----------------------------- 
    
    relf  = f/f0; %Relative reduction in objective f
    merit = ps;


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

   if (mod(itn,pars.itn_print) == 0 || itn ==1)&&(pars.itn_print>0)
            fprintf('%5.1i %3.1i %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e\n',...
                     itn,k,relres,relf,merit,mu,nP/nP0,nD/nD0,nG/nG0,t,min(x),min(z),infeas,rel_infeas); 
   end

   
    %---------------------- 
    %Form the gradient of psi
    %----------------------
    sc_l = 2*rho/(nP^2 + nD^2 + nG^2);
    grad = sc_l*2*[A'*p+g*c;A*d-g*b;d] - [1./x;zeros(m,1);1./z]; 

    
    %Save the present step    
    xs = x;
    ys = y;
    zs = z;
    
    %-----------------------------------
    %Backtracking linesearch
    %----------------------------------

    for k = 1:pars.max_bkt
        max_step = [xs./grad(1:n);zs./grad(n+m+1:2*n+m)];
        pos_ix   = find(max_step >=0);
        max_step = min(max_step(pos_ix));
        max_step = max_step*.99;
        
        t        = min(t,max_step);

        %-------------------------------
        %Take the gradient step
        %-------------------------------
        x  = xs - t*grad(1:n);
        y  = ys - t*grad(n+1:n+m);
        z  = zs - t*grad(n+m+1:2*n+m); 

        nGradS = norm(grad)^2;

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
       
        %Calculate the achieved decrease
        fn  = nP^2+nD^2+nG^2;
        psn = 2*rho*log(fn) - sum(log(x)) - sum(log(z));
        if(~isreal(psn))
            'notreal'
        end

 
        %Evaluate the backtracking inequality 
        if psn < ps - t*0.5*nGradS
             break;
        else %Backtracking step
            t = t*pars.beta_bkt;
        end


    end % end of linesearch

    f = fn;
    ps= psn;

    if(stop_in_error)
      results.itn = itn;
      break;
    end
    %---------------------------------
    %Nesterov/ P.Tseng Acceleration
    %--------------------------------- 
    if pars.accel > 0  %Acceleration step
        theta   = (itn-1)/(itn+2); %Extrapolation coeficient
        ex_dirx = [x-xp];
        ex_diry = [y-yp];
        ex_dirz = [z-zp];

        grad    = [ex_dirx;ex_diry;ex_dirz];

        %Save x
        xp      = x;
        zp      = z;
        yp      = y;
        resp    = res;
 

        %Save the present step    
        xs = x;
        ys = y;
        zs = z;
 
        %-----------------------------------
        %Backtracking linesearch
        %----------------------------------

        for k = 1:pars.max_bkt
            max_step = [xs./grad(1:n);zs./grad(n+m+1:2*n+m)];
            pos_ix   = find(max_step >=0);
            max_step = min(max_step(pos_ix));
            max_step = max_step*.99;
            
            t_ep        = min(t_ep,max_step);

            %-------------------------------
            %Take the gradient step
            %-------------------------------
            x  = xs - t_ep*grad(1:n);
            y  = ys - t_ep*grad(n+1:n+m);
            z  = zs - t_ep*grad(n+m+1:2*n+m); 

            nGradS = norm(grad)^2;

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
           
            %Calculate the achieved decrease
            fn  = nP^2+nD^2+nG^2;
            psn = 2*rho*log(fn) - sum(log(x)) - sum(log(z));
            if(~isreal(psn))
                'notreal'
            end

     
            %Evaluate the backtracking inequality 
            if psn < ps - t*0.5*nGradS
                 break;
            else %Backtracking step
                t_ep = t_ep*pars.beta_bkt;
            end


        end % end of linesearch

        ps = psn;

    end

 
   %Prepare the next loop
   if(isfield(pars,'t0m')) %Extend the initial step size a little
        t = t*pars.t0m;
        t_ep = t_ep*pars.t0m;

   end
          
end

%---------------------------------------------------------
%End of the iteration
%---------------------------------------------------------
results.itn = itn;
