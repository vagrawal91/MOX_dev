clear all;
close all;

load data/genes.mat
X = X-mean(X);
Y = zscore(Y);

m = size(X, 1);
p = size(X, 2);
q = size(Y, 2);
r = min(p,q);


l = 3;
[P, D, Q, muX, muY, E, Fmax, A, B, W] = mox(X, Y, r, l);
% Sort genes accoring to a Euclidean norm
nP = sqrt(sum(P.^2, 2));
[~,idx] = sort(-nP);
P = P(idx,:);
X = X(:,idx);

fig1 = figure(1);
% Create a tiled layout instead of subplots
t = tiledlayout(3, 6, 'TileSpacing', 'Compact', 'Padding', 'Compact');
ax1 = nexttile(t, [1 3]); % Spanning 3 tiles
    bar(1e3*(D*Q')');
    xticks(1:size(Q, 1));
	xticklabels(drugNames);
	ylabel('{\bfDQ}^T [10^{-3}]');
    subplotlabel('a', 'northwest');

ax2 = nexttile(t, [3 3]); % Spanning 3 tiles

% Retrieve MATLAB's default color order
defaultColors = get(groot, 'defaultAxesColorOrder');

% Define clusters with indices, using default colors and specified markers and labels
clusters = {
  struct('indices', [7, 10], 'color', defaultColors(1, :), 'marker', 'o', 'label', 'Microtubule Inhibitors', 'valign', 'top', 'halign', 'left'),   % Paclitaxel, Vinblastine
    struct('indices', [5, 6, 9], 'color', defaultColors(2, :), 'marker', 's', 'label', 'Topoisomerase Inhibitors', 'valign', 'bottom', 'halign', 'left'), % Etoposide, Irinotecan, Topotecan
    struct('indices', [3, 2], 'color', defaultColors(3, :), 'marker', '^', 'label', 'DNA Cross-Linkers', 'valign', 'top', 'halign', 'left'),    % Cisplatin, Carboplatin
    struct('indices', [4], 'color', defaultColors(4, :), 'marker', 'v', 'label', 'DNA Intercalating Agent', 'valign', 'bottom', 'halign', 'right'),       % Doxorubicin
    struct('indices', [1], 'color', defaultColors(5, :), 'marker', '>', 'label', 'DNA Cleaving Agent', 'valign', 'top', 'halign', 'right'),       % Bleomycin
    struct('indices', [8], 'color', defaultColors(6, :), 'marker', 'd', 'label', 'Hormonal Agent', 'valign', 'top', 'halign', 'right')             % Tamoxifen
};

% Initialize legend entries
legendEntries = {};

% Plot each cluster with its specific color, marker, and add to legend
hold on;
hs = [];
for c = 1:length(clusters)
    % Plot each point with specified color and marker
	  h=plot3(Q(clusters{c}.indices, 1), Q(clusters{c}.indices, 2), Q(clusters{c}.indices, 3), ...
          clusters{c}.marker, 'MarkerSize', 8, 'MarkerFaceColor', clusters{c}.color, ...
          'MarkerEdgeColor', 'none', 'LineStyle', 'none');
    legendEntries{end+1} = clusters{c}.label; % Add label to legend entries
    hs = [hs h];

    % Draw dotted lines and black dots at the intersection with the plane
    for i = clusters{c}.indices
        plot3([Q(i, 1), Q(i, 1)], [Q(i, 2), Q(i, 2)], [Q(i, 3), 0], ...
              'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        plot3(Q(i, 1), Q(i, 2), 0, 'ko', 'MarkerSize', 1.5, 'MarkerFaceColor', 'k');
        text(Q(i, 1), Q(i, 2), Q(i, 3), ['  ' drugNames{i}], 'VerticalAlignment', clusters{c}.valign, 'HorizontalAlignment', clusters{c}.halign, 'fontsize', 10);
    end
end

% Plot a semi-transparent reference plane at z = 0
mrg = 0.15;
xLimits = [(min(Q(:,1))-mrg) (max(Q(:,1))+mrg)];
yLimits = [(min(Q(:,2))-mrg) (max(Q(:,2))+mrg)];
zLimits = [(min(Q(:,3))-mrg) (max(Q(:,3))+mrg)];
[X1, Y1] = meshgrid(linspace(xLimits(1), xLimits(2), 10), linspace(yLimits(1), yLimits(2), 10));
Z1 = zeros(size(X1));
plane = surf(X1, Y1, Z1, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);

% Set plot title and labels
xlabel('{\bfQ}_{(:,1)}');
ylabel('{\bfQ}_{(:,2)}');
zlabel('{\bfQ}_{(:,3)}');
box on;
view(-50, 22);
axis([xLimits yLimits zLimits]);
subplotlabel('d', 'northeast');

% Add legend with cluster labels
legend(hs, legendEntries, 'Location', 'northwest', 'fontsize', 8);
hold off;

	
% Plot (prediction)
ax3 = nexttile(t, [2 2]);
    mkr = {'o';'s';'d';'+';'x';'*';'^';'<';'v';'>';'o';'s';'d';'+';'x';'*';'^';'<';'v';'>'};
    Yhat = X*P*D*Q';
    h = plot(Y, Yhat, '+');
    set(h, {'Marker'}, mkr(1:size(Yhat,2)));
    hold on;
    plot(10*[-1 1], 10*[-1 1], 'k-');
    xline(0, ':k');
    yline(0, ':k');
    hold off;
    axis([-5 5 -6 4]);
    set(ax3, 'layer', 'top', 'box', 'on');
    xticks(-10:10);
    xlabel('{\bfY}');
	ylabel('{\bfY}_{pred}');
	lgd2 = legend(drugNames, 'location', 'south', 'numcolumns', 2);
    subplotlabel('b', 'northwest');

	

    % Fifth plot (residual histogram)
    ax4 = nexttile(t);
    res = Y - Yhat;
    res = res(:);
    bins = round(log2(length(res)) + 1);
    if mod(bins, 2) == 0
      bins = bins + 1;
    end
    x_lim = max(abs(res)*(1+0.5/bins));
    histogram(res, linspace(-x_lim, x_lim, bins+1), 'Normalization', 'pdf');
    hold on;
    pd = fitdist(res, 'Normal');
    x_values = linspace(-x_lim, x_lim, 200);
    y_values = pdf(pd, x_values);
    plot(x_values, y_values, '-', 'LineWidth', 2);
    xlim(x_lim*[-1 1]);
    xlabel('{\bfY}_{pred} - {\bfY}');
    ylabel('Probability Density');
    hold off;
    subplotlabel('c', 'northwest');
	
    % Sixth plot: Information text tile
    ax5 = nexttile(t);
    % Compute R-squared for the fit
    SS_res = sum((Y - Yhat).^2, 'all');
    SS_tot = sum((Y - mean(Y)).^2, 'all');
    R2 = 1 - SS_res / SS_tot;

	
    % Display dataset name, number of latent variables, and R^2
    text(0.01, 0.7, 'Dataset: Gene', 'FontSize', 8, 'FontWeight', 'bold');
    text(0.01, 0.5, ['Latent variables: ', num2str(l)], 'FontSize', 8);
    text(0.01, 0.3, ['{\itR}^2: ', sprintf('%.4f', R2)], 'FontSize', 8);
    set(ax5, 'Visible', 'off'); % Hide axes for text-only tile

fsize = [29 13]; % cm
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 fsize], 'PaperSize', fsize);
set(gcf, 'Units', 'centimeters', 'Position', [0 0 fsize]);
exportgraphics(gcf, 'img/mox_gene.pdf', 'ContentType', 'vector');

