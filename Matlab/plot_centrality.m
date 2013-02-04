%Plots the centrality x'z/n of function f1, f2 extracted from the lpnetlib results
%Requires the variable prob_index defined, and the cell arrays resultsf1 and results f2

meritf1 = resultsf1{prob_index}{11}{1};
meritf2 = resultsf2{prob_index}{11}{1};
figure
hold on
plot(meritf1,'bd');
plot(meritf2,'rd');
legend({'f1','f2'});
hold off



