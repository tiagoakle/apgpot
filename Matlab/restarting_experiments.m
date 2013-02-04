clear all;
clc;
%December 16 2012
%iterates over the lpnetlib problems in standard form and
%uses the objectives f1 and f2 to solve. For each function and each problem
%it records problem name, problem size , iteration count, norm of primal residual
%norm of dual residual, centrality and gap.

%Load the indices of the standard form 
%problems in lpnetlib
resultsf1 = {{'Problem Name','Objective function','m','n','Iterations','Primal Res Norm','Dual Res Norm','Gap','Centrality','exit flag'}};
resultsf2 = {{'Problem Name','Objective function','m','n','Iterations','Primal Res Norm','Dual Res Norm','Gap','Centrality','exit flag'}};
load 'standard_form_indices';
problem_count = length(st_ix);
for problem_index = 1:problem_count
  problem_uf_ix = st_ix(problem_index);
  P = UFget(problem_uf_ix);

  prob_name = [P.name];
  prob_name(find(prob_name=='/'))=' ';

  A = P.A;
  %Problem parameters
  m = size(A,1);
  n = size(A,2);

  b = P.b;
  c = P.aux.c;
  fprintf('Loaded problem %s, m:%i, n:%i\n',P.name,m,n);


  pars.rho = n+sqrt(n);
  pars.accel = 1;
  pars.max_iter = min(m*n,30000);
  pars.tolnorm = 2;
  pars.itn_print = 30000;
  pars.restart   = 5000;

  %Run prox_pot_toy to compare
  prox_pot_toy
  log_cell = {log_mu,...
              log_p_r,...
              log_d_r,...   
              log_g_r   ,...
              log_merit ,...
              log_f     ,...
              log_step  ,...
              log_infeas,...
              log_obj};
  res_cell = {prob_name,'1',m,n,itn,nP,nD,nG,mu,Exit_flag,log_cell};
  resultsf1{problem_index+1} = res_cell;

  prox_pot_toy_skew
  log_cell = {log_mu,...
              log_p_r,...
              log_d_r,...   
              log_g_r   ,...
              log_merit ,...
              log_f     ,...
              log_step  ,...
              log_infeas,...
              log_obj};
  res_cell = {prob_name,'2',m,n,itn,nP,nD,nG,mu,Exit_flag,log_cell};
  resultsf2{problem_index+1} = res_cell;
  save 
end
