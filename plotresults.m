% This is a MATLAB/Octave script.
% Run this script using command "octave -q -f --eval plotresults" after 
% obtaining results.m by running hw2_6.

% Run the results script
results;
N = size(digits_cgs,1)+1;

% Plot errors and print into file
figure(1);
range = 2:N;
plot(range, digits_cgs, 'r-o', range, digits_mgs, 'b--x', range, digits_dgs, 'k-.');
title('Digits of accuracy');
xlabel('dimension');
ylabel('digits of accuracy');
legend('CGS','MGS','DCG','Location','NorthEast');
print -deps errors.eps

% Plot execution times and print into file
figure(2);
plot(range, times_cgs, 'r-o', range, times_mgs, 'b--x', range, times_dgs, 'k-.');
title('Execution times');
xlabel('dimension');
ylabel('wallclock time in seconds');

legend('CGS','MGS','DCGS','Location','NorthWest');
print -deps times.eps
