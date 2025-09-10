clear all;
close all;

% Plot the average cross-validated MSE as a function of the number of components
fig1 = figure(1);
t = tiledlayout(1, 5, 'TileSpacing', 'Compact', 'Padding', 'Compact');
lblfontsize = 12;

load result_simulated_dataset.mat
ax1 = nexttile(t, [1 2]); % Spanning 3 tiles
hold on;
errorbar(r*[-1 0.5 2], MSEols_mean*[1 1 1], MSEols_std*[1 1 1], '--', 'DisplayName', 'OLS');
errorbar(1:r, MSEcca_mean, MSEcca_std, '-o', 'DisplayName', 'CCA');
errorbar(1:r, MSEpls_mean, MSEpls_std, '-s', 'DisplayName', 'PLS');
errorbar(1:r, MSEmox_ll_mean, MSEmox_ll_std, '-^', 'DisplayName', 'MOX_{{\ithh}}');
errorbar(1:r, MSEmox_rl_mean, MSEmox_rl_std, '-d', 'DisplayName', 'MOX_{{\itrh}}');
hold off;
xlabel('\ith');
ylabel('Cross-validated MSE');
title('Simulated Dataset')
legend('Location', 'northeast');
xlim([0 r]);
ylim([0 2]);
subplotlabel('a', 'northwest');
ax = gca;
xticks(0:r);
%set(ax, 'XMinorTick', 'on');
%ax.XAxis.MinorTickValues = 0:r;
set(ax, 'layer', 'top', 'box', 'on');

load result_corn_dataset.mat
ax1 = nexttile(t, [1 1]);
hold on;
plot(r*[-1 0.5 2], MSEols_mean*[1 1 1], '--', 'DisplayName', 'OLS');
plot(1:r, MSEcca_mean, '-o', 'DisplayName', 'CCA');
plot(1:r, MSEpls_mean, '-s', 'DisplayName', 'PLS');
plot(1:r, MSEmox_ll_mean, '-^', 'DisplayName', 'MOX_{{\ithh}}');
plot(1:r, MSEmox_rl_mean, '-d', 'DisplayName', 'MOX_{{\itrh}}');
hold off;
xlabel('\ith');
ylabel('Cross-validated MSE');
title('Corn Dataset')
legend('Location', 'northeast');
xlim([0 r]);
ylim([0 3]);
subplotlabel('b', 'northwest');
ax = gca;
xticks(0:r);
set(ax, 'layer', 'top', 'box', 'on');

load result_pulp_dataset.mat
ax1 = nexttile(t, [1 1]);
hold on;
plot(r*[-1 0.5 2], MSEols_mean*[1 1 1], '--', 'DisplayName', 'OLS');
plot(1:r, MSEcca_mean, '-o', 'DisplayName', 'CCA');
plot(1:r, MSEpls_mean, '-s', 'DisplayName', 'PLS');
plot(1:r, MSEmox_ll_mean, '-^', 'DisplayName', 'MOX_{{\ithh}}');
plot(1:r, MSEmox_rl_mean, '-d', 'DisplayName', 'MOX_{{\itrh}}');
hold off;
xlabel('\ith');
ylabel('Cross-validated MSE');
title('Pulp Dataset')
legend('Location', 'northeast');
xlim([0 r]);
ylim([0 1.6]);
subplotlabel('c', 'northwest');
ax = gca;
xticks(0:2:r);
set(ax, 'XMinorTick', 'on');
ax.XAxis.MinorTickValues = 0:r;
set(ax, 'layer', 'top', 'box', 'on');


load result_gene_dataset.mat
ax1 = nexttile(t, [1 1]);
hold on;
plot(1:r, MSEpls_mean, '-s', 'DisplayName', 'PLS');
plot(1:r, MSEmox_rl_mean, '-d', 'DisplayName', 'MOX_{{\itrh}}');
hold off;
xlabel('\ith');
ylabel('Cross-validated MSE');
title('Gene Dataset')
legend('Location', 'northeast');
xlim([0 r]);
ylim([0 1.6]);
subplotlabel('d', 'northwest');
ax = gca;
xticks(0:2:r);
set(ax, 'XMinorTick', 'on');
ax.XAxis.MinorTickValues = 0:r;
set(ax, 'layer', 'top', 'box', 'on');




fsize = [28 7]; % cm
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 fsize], 'PaperSize', fsize);
set(gcf, 'Units', 'centimeters', 'Position', [0 0 fsize]);
%exportgraphics(gcf, 'img/dimensionality.pdf', 'ContentType', 'vector');
