%January 17 2013
%This script generates a random problem and 
%executes prox_pot_toy with as many different values of rho as the cell array rhos contains
%It then plots the value of f for each iteration 
clear all;
close all;
clc;
addpath ./apgpot
%Generate a random problem
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

rhos = { n+sqrt(n),'r','n + sqrt(n)';...
         2*n      ,'g' ,'2n';...    
         4*n      ,'b' ,'4n';...
         m*n      ,'k' ,'mn';...
         n*n      ,'y' ,'n^2';...
         n*n*n    ,'c' ,'n^3';...
         n*n*n*n  ,'m' ,'n^4'};


pars.accel = 1;
h = figure;
hold on;
for r = 1:size(rhos,1)
  pars.rho   = rhos{r,1};
  %--------------------------------------------------------
  %Run prox_pot_toy with the potential reduction projection
  %--------------------------------------------------------
  pars.proj = 'pot';
  prox_pot_toy;
  %Save the results
  f_pot = results.log_f(1:results.itn);
  plot(log10(f_pot),rhos{r,2});
end
legend(rhos(:,3));
saveas(h,'../Writeup/Images/potvsrho.png');
hold off

