clear all;
close all;

%load result_simulated_dataset.mat
%load result_pulp_dataset.mat
%load result_corn_dataset.mat

% Plot the average cross-validated MSE as a function of the number of components
fig1 = figure(1);
%subplot(1, 2, 1);
hold on;
errorbar(r*[-1 0.5 2], MSEols_mean*[1 1 1], MSEols_std*[1 1 1], '--', 'DisplayName', 'OLS');
errorbar(1:r, MSEcca_mean, MSEcca_std, '-o', 'DisplayName', 'CCA');
errorbar(1:r, MSEpls_mean, MSEpls_std, '-s', 'DisplayName', 'PLS');
errorbar(1:r, MSEmox_ll_mean, MSEmox_ll_std, '-^', 'DisplayName', 'MOX_{ℓℓ}');
errorbar(1:r, MSEmox_rl_mean, MSEmox_rl_std, '-d', 'DisplayName', 'MOX_{{\itr}ℓ}');
hold off;
xlabel('Number of Components (ℓ)');
ylabel('Cross-validated MSE');
legend('Location', 'northeast');
xlim([0 r]);
%ylim([0.0 * min([MSEmox_mean, MSEpls_mean, MSEols_mean]), ...
%      2.0 * max([MSEmox_mean, MSEpls_mean, MSEols_mean])]);
ylim([0 2]);
ax = gca;
set(ax, 'layer', 'top', 'box', 'on');
currentPosition = fig1.Position; 
fig1.Position = [currentPosition(1), currentPosition(2), 0.9*currentPosition(3), 0.9*currentPosition(4)];
		      
% Set up printing options
set(fig1, 'Color', 'white');   % Ensure the figure background is white (optional)
set(fig1, 'Renderer', 'painters');
print(fig1, '-depsc', 'img/benchmark_simulated.eps');
