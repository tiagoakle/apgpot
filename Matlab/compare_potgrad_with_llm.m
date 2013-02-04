%Jan 30 2013
%The objective of this experiment is to compare the behavior of the protential gradient
%proximal potential approach with the behavior of the LLM approach.

addpath ../Matlab/TfocsTests %Add the path to the tfocs_llm_compare code
clear all;
close all;

%Define problem parameters
m=5;
n=15;
density = 5;


  %Generate the problem instance
  A = sprand(m,n,density/100);
  A = randn(m,n);
  x0    = rand(n,1);
  b     = A*x0; %Find a rhs which makes the problem feasible 
  y0    = randn(m,1);
  z0    = rand(n,1);
  c     = A'*y0 + z0; %Generate a feasible dual

  fprintf(' %4s %4s %6s %6s %6s %6s %6s \n','n','m','spa','t_itn','t_tim','p_itn','p_tim');

  %---------Call TFOCS LLM----------------------------------------------------------
  clear opts
  opts = struct;
  opts.alg = 'LLM';
  opts.maxIts = 600;

  E = [A      ,zeros(m,m),zeros(m,n) ;...
      zeros(n),A'        ,eye(n)     ;...
      c'      ,-b'       ,zeros(1,n)];
  q = [-b;-c;0];
  tic
  [sol,res] = tfocs(smooth_quad,{E,q},proj_restrictedRplus(n,m),ones(2*n+m,1),opts);
  tfocs_time = toc;
  tfocs_res = res;
  tfocs_niter = res.niter;
  legends = {'LLM Tfocs'};
  fprintf(' %4i %4i %4i %6i %6f \n',n,m,density,tfocs_niter,tfocs_time);
  %------End of tfocs call ---------------------------------------------------------
    
  %print 
  close all;
  semilogy(res.f,'c');
  hold on;

%-----------------Call prox pot toy llm-------------------------------------------
pars.proj = 'ind'; %The projection must be the indicator
pars.max_iter = 600;
pars.accel = 1;
prox_pot_toy;
%Save the results
f_ind = results.log_f(1:results.itn);
semilogy(f_ind,'b');
legends = {legends{:},'LLM mine'};
%-----------------End of call to prox_pot_toy------------------------------------

%-----------------Call prox pot toy llm-------------------------------------------
pars.proj = 'pot'; %The projection must be the potential
pars.max_iter = 600;
pars.accel = 1;
pars.rho = n+sqrt(n);
prox_pot_toy;
%Save the results
f_pot = results.log_f(1:results.itn);
semilogy(f_pot,'r--');
legends = {legends{:},'prox_pot'};
%-----------------End of call to prox_pot_toy------------------------------------


 %----------------Call potgrad_prox_test with different values of rho-------------
 rhos = {n+1,'n+1','k';
 n+sqrt(n),'n+sqrt(n)','r';
 n+m,'n+m','m';
 2*n,'2n','g';
 n*n,'n^2','y'};
 
 for rhi_ix = 1:size(rhos,1)
   pars.rho = rhos{rhi_ix,1}; %Set the parameter
   pars.max_iter = 600;       %Set the iteration bound
   pars.accel = 1;            %Set the acceleration to true
   tic                        %Time the process
   potgrad_prox_test;
   potgrad_time = toc;  
   potgrad_niter = results.itn; 
   fprintf(' %6i %6f %6f \n',potgrad_niter,potgrad_time,rho);
   semilogy(results.log_f(1:results.itn),rhos{rhi_ix,3})
   legends = {legends{:},rhos{rhi_ix,2}};
 end
 %----------------End of call to potgrad prox test ------------------------------

  legend(legends);
  hold off;
hold off;
