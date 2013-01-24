%January 16 2013
%This script generates a random problem and 
%executes prox_pot_toy twice, the first time with pars.proj = 'pot' and the second with pars.proj = 'ind'
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


%%--------------------------------------------------------
%%Run prox_pot_toy with the potential gradient
%%--------------------------------------------------------
%pars.proj = 'pot';
%pars.func = 'psi';
%pars.max_iter =20000;
%prox_pot_toy;
%%Save the results
%f_ind = results.log_f(1:results.itn);


pars.accel = 1;
pars.rho   = n+sqrt(n);
%--------------------------------------------------------
%Run prox_pot_toy with the potential reduction projection
%--------------------------------------------------------
pars.proj = 'pot';
prox_pot_toy;
%Save the results
f_pot = results.log_f(1:results.itn);

%--------------------------------------------------------
%Run prox_pot_toy with the indicator projection
%--------------------------------------------------------
pars.proj = 'ind';
pars.max_iter =20000;
prox_pot_toy;
%Save the results
f_ind = results.log_f(1:results.itn);

%--------------------------------------------------------
%Run lpapgpot to compare
%--------------------------------------------------------
clear pars 
pars = struct;
pars.maxit = 20000;
pars.tol   = 1.e-6;
[x,f,R]=lpapgpot(A,b,c,pars);
%Save the results
f_apg = R.log_f(1:R.iter);


%-------------------------------------------------------
%Plot the results
%-------------------------------------------------------
h = figure;
hold on
plot(log10(f_pot),'r');
plot(log10(f_ind),'b');
legend('Potential projection','Indicator projection');
saveas(h,'../Writeup/Images/potvsprox.png');
hold off

%-------------------------------------------------------
%Plot the results
%-------------------------------------------------------
h = figure;
hold on
plot(log10(f_pot),'r');
plot(log10(f_ind),'b');
plot(log10(f_apg),'g');
legend('Potential projection','Indicator projection','APGPOT');
saveas(h,'../Writeup/Images/potvsproxvsapgpot.png');
hold off
