%Jan 25 2012:
%Compare the rate of convergence for the models
%Pr: min c^Tx + 1/(2mu)||Ax-b||_2^2 st 0<=x
%Dr: min -b^Ty + 1/(2mu)||A^Ty+z-c||_2^2 st 0<=z
%PDG: min (c^Tx -b^Ty) + ||Ax-b||_2^2 + ||A^Ty+z-c||_2^2  st 0<= x,z

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


%-------Solve primal dual gap------------------------
Epdg = [A      ,zeros(m,m),zeros(m,n) ;...
    zeros(n),A'        ,eye(n)     ;...
    c'      ,-b'       ,zeros(1,n)];

opts = struct;
opts.saveHist = true;
opts.maxIts   = 1000;
opts.alg      = 'AT'
opts.errFcn   = @(f,x)(c'*x(1:n)); %function to store the objective value of the lp
[sol_pdg,res_pdg,out_pars] = tfocs(smooth_quad,{Epdg},proj_restrictedRplus(n,m),ones(2*n+m,1),opts);
%------------End of primal dual gap ------------------

%smu = sqrt(mu) 
mu = 1.e-10;
smu = sqrt(mu);
%------------Solve with the regularized primal--------
smootF = {smooth_linear(c),smooth_quad};
affineF={1/smu*A,-1/smu*b;c,0}
opts.errFcn   = @(f,x)(c'*x(1:n)); %function to store the objective value of the lp
[sol_preg,res_preg,out_pars] = tfocs(smooth_quad,affineF,proj_Rplus,[],opts);
%------------End of regularized primal----------------

%-----------Solve with regularized dual---------------
affineF={[1/smu*A',1/smu*eye(n)],-c;[-b;zeros(n,1)],0};
[sol_dreg,res_dreg,out_pars] = tfocs(smooth_quad,affineF,proj_restrictedRplusDual(n,m),[],opts);
%-----------End of regularized dual

hold on
plot(res_pdg.err,'bd');
plot(res_preg.err,'r');
plot(res_dreg.err,'g');
hold off
