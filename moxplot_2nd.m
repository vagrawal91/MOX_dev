function [axs, lgds] = moxplot_2nd(P, D, Q, X, Y, datasetName, Pxlabels, DQxlabels, Pylabel, DQylabel, plotAsSpectra, lambda)
% MOXPLOT Plots MOX regression results as bar charts or line plots for spectra.
% 
% This function creates either bar plots or line plots (for spectra) for the
% predictor loadings (P), impact (D), and response loadings (Q) in a tiled layout.
% The function allows for customizable x-axis and y-axis labels, but uses default labels if none are provided.
%
% Additional feature: Adds scientific paper-style letter labels (a, b, c, etc.) to each subplot.
%
% USAGE:
%   [axs, lgds] = moxplot(P, D, Q, X, Y, datasetName)
%   [axs, lgds] = moxplot(P, D, Q, X, Y, datasetName, Pxlabels, DQxlabels, Pylabel, DQylabel, plotAsSpectra, lambda)
%
% INPUTS:
%   P             - Predictor loadings matrix 
%   D             - Diagonal impact matrix
%   Q             - Response loadings matrix
%   X             - Predictor matrix
%   Y             - Response matrix
%   datasetName   - Name of the dataset to be displayed
% 
% Optional Inputs:
%   Pxlabels      - Cell array of x-axis labels for predictor loadings (default: integer strings)
%   DQxlabels      - Cell array of x-axis labels for response loadings (default: integer strings)
%   Pylabel       - String for y-axis label of the predictor loadings plot (default: '{\bfP}')
%   DQylabel       - String for y-axis label of the response loadings plot (default: '{\bfQ}')
%   plotAsSpectra - If true, plots P as spectra using lines; default is false (bar chart)
%   lambda        - Wavelengths for spectral plot; required if plotAsSpectra is true
%
% OUTPUT:
%   axs           - Vector of axis handles [ax1, ax2, ax3, ax4].
%   lgds          - Handle for the legends [lgd1, lgd2].
%

    nLatent = size(P, 2);
    lblfontsize = 12;
  
    % Default labels if not provided
    if nargin < 7 || isempty(Pxlabels)
        Pxlabels = string(1:size(P, 1)); % Default Pxlabels as integer strings
    end
    if nargin < 8 || isempty(DQxlabels)
        DQxlabels = string(1:size(Q, 1)); % Default Qxlabels as integer strings
    end
    if nargin < 9 || isempty(Pylabel)
        Pylabel = '{\bfP}'; % Default Pylabel
    end
    if nargin < 10 || isempty(DQylabel)
        DQylabel = '{\bfDQ}^T'; % Default DQylabel
    end
    if nargin < 11 || isempty(plotAsSpectra)
        plotAsSpectra = false; % Default is to plot as bar chart
    end

    % Check for lambda if plotAsSpectra is true
    if plotAsSpectra && (nargin < 12 || isempty(lambda))
        error('Wavelengths (lambda) must be provided for spectral plots.');
    end

    % Create a tiled layout instead of subplots
    t = tiledlayout(4, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    % First plot (P)
    ax1 = nexttile(t, [1 3]); % Spanning 3 tiles
    if plotAsSpectra
        % Plot as spectra (lines) with wavelengths on x-axis
        plot(lambda, cumsum(cumsum(P, 2), 2), 'LineWidth', 2);
        xlabel('Wavelength (\lambda)');
	xlim([floor(min(lambda)/100)*100 ceil(max(lambda)/100)*100]);
    else
        % Plot as bar chart
        bar(P);
        xticks(1:size(P, 1));
        xticklabels(Pxlabels);
    end
    xlabel('Predictor');
    ylabel(Pylabel);
    
    % Dynamically create the legend based on the number of latent variables (columns in P)
    latentLabels = arrayfun(@(i) sprintf('latent %d', i), 1:nLatent, 'UniformOutput', false);
    
    % Add a horizontal legend spanning the width of the plot
    lgd1 = legend(latentLabels, 'Orientation', 'horizontal');
    lgd1.Layout.Tile = 'north';  % Position the legend above the first tile
    subplotlabel('a', 'northwest');

    % Plot (DQ')
    ax2 = nexttile(t, [1 3]); % Spanning 2 tiles
    bar((D*Q')');
    xticks(1:size(Q, 1));
    xticklabels(DQxlabels);
    ylabel(DQylabel);
    subplotlabel('b', 'northwest');

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
    sc = ceil(max(abs([Y(:);Yhat(:)])));
    axis(sc*[-1 1 -1 1]);
    set(ax3, 'layer', 'top', 'box', 'on');
    xlabel('{\bfY}');
    ylabel('{\bfY}_{pred}');
    lgd2 = legend(DQxlabels, 'location', 'southeast');
    subplotlabel('c', 'northwest');

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
    subplotlabel('d', 'northwest');
    
    % Sixth plot: Information text tile
    ax5 = nexttile(t);
    % Compute R-squared for the fit
    SS_res = sum((Y - Yhat).^2, 'all');
    SS_tot = sum((Y - mean(Y)).^2, 'all');
    R2 = 1 - SS_res / SS_tot;
    
    % Display dataset name, number of latent variables, and R^2
    text(0.01, 0.7, ['Dataset: ', datasetName], 'FontSize', 8, 'FontWeight', 'bold');
    text(0.01, 0.5, ['Latent variables: ', num2str(nLatent)], 'FontSize', 8);
    text(0.01, 0.3, ['{\itR}^2: ', sprintf('%.4f', R2)], 'FontSize', 8);
    set(ax5, 'Visible', 'off'); % Hide axes for text-only tile

    % Return the axis handles and legend handles
    axs = [ax1, ax2, ax3, ax4, ax5];
    lgds = [lgd1, lgd2];
end


