%Plots the primal residual of function f1, f2 extracted from the lpnetlib results
%Requires the variable prob_index defined, and the cell arrays resultsf1 and results f2

dresf1 = resultsf1{prob_index}{11}{3};
dresf2 = resultsf2{prob_index}{11}{3};
figure
hold on
plot(dresf1,'bd');
plot(dresf2,'rd');
legend({'f1','f2'});
hold off
