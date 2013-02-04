%Jan 30 2013
%Generates problems comparable to those
%reported by llm and proceeds to solve them with TFOCS using the file tfocs_llm.m
%For each problem it saves the iteration count, and the full results file.
clear all
%Reset the seed
seed = 0;
randn('seed',seed);
rand('seed',seed);

fprintf(' %4s %4s %4s %6s %6s\n','n','m','dens','niter','time');
fprintf('---------------------------------- \n');
results = {};
%Main loop
diary('results.txt');
for problem_count =1:21
  %Generate the problem parameters
  n=0;
  m=0;
  density=0;

  if(problem_count<=9)
    n=1000;
    m=floor((problem_count-1)/3)+1;
    if(m==1);m=100;end;
    if(m==2);m=500;end;
    if(m==3);m=900;end;
  elseif(problem_count<=18)
    n=5000;
    m=floor((problem_count-1)/3)+1;
    if(m==4);m=500;end;
    if(m==5);m=2500;end;
    if(m==6);m=4500;end;
  else
    n=10000;
     m=problem_count;
    if(m==19);m=1000;end;
    if(m==20);m=5000;end;
    if(m==21);m=9000;end;
  end

  if(problem_count<=18)
    density = 5*(mod(problem_count-1,3));
    if(density==0); density=1; end;
  else
    density = 1;
  end

  %Generate the problem instance
  A = sprand(m,n,density/100);

  x0    = rand(n,1);
  b     = A*x0; %Find a rhs which makes the problem feasible
  
  y0    = randn(m,1);
  z0    = rand(n,1);
  c     = A'*y0 + z0; %Generate a feasible dual
  %End of instance generation
  opts = struct;
  opts.printEvery = 0; %supress output
  tic;
  tfocs_llm_compare
  time =toc;
  res.time = toc;
  %Save the results
  niter = res.niter;
  fprintf(' %4i %4i %4i %6i %6i \n',n,m,density,niter,time);
  results = {results{:},res};
  save 'results' 'results' %Save the results in results mat
  
end
    
diary off;


