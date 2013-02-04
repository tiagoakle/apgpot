%Feb 1 2013 based of compare_potgrad_with_llm
%The objective of this experiment is to compare the behavior of the protential gradient
%proximal potential approach with the behavior of the LLM approach, and the proximal potential approach.

addpath ../Matlab/TfocsTests %Add the path to the tfocs_llm_compare code
clear all;
close all;

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


%common parameters 
max_iter  = 600;

%--------------------------------------------------------
%Run prox_pot_toy with the potential reduction projection
%--------------------------------------------------------
pars.accel = 1;
pars.rho   = n+sqrt(n);
pars.proj = 'pot';
pars.max_iter =max_iter;

prox_pot_toy;
%Save the results
f_pot = results.log_f(1:results.itn);

%--------------------------------------------------------
%Run prox_pot_toy with the indicator projection
%--------------------------------------------------------
pars.accel = 1;
pars.rho   = n+sqrt(n);
pars.proj = 'ind';
pars.max_iter =max_iter;

prox_pot_toy;
%Save the results
f_ind = results.log_f(1:results.itn);

%-------------------------------------------------------
%Run potgrad test
%-------------------------------------------------------
pars.accel = 1;
pars.rho   = n+sqrt(n);
pars.max_iter =max_iter;
potgrad_prox_test;
f_potgrad = results.log_f(1:results.itn);

%------------------------------------------------------
%Print results
%------------------------------------------------------
h = semilogy(f_pot,'bd');
hold on;
semilogy(f_ind,'r-.');
semilogy(f_potgrad,'g-');
legend('Prox Pot','LLM','Prox PotGrad');
hold off;
saveas(h,'../Writeup/Images/llm_potgrad_proxpot.png');
