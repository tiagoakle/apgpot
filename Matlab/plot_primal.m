%Plots the primal residual of function f1, f2 extracted from the lpnetlib results
%Requires the variable prob_index defined, and the cell arrays resultsf1 and results f2

presf1 = resultsf1{prob_index}{11}{2};
presf2 = resultsf2{prob_index}{11}{2};
figure
hold on
plot(presf1,'bd');
plot(presf2,'rd');
legend({'f1','f2'});
hold off
