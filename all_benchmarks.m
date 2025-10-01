clear; clc; close all;
clearvars -except filename m p q w_amp n_random_seeds n_repetitions CV d_XY ...
    eps_x eps_y d_X d_Y d idx_X idx_XY idx_Y iex iey im iq ip ...
    m_vec p_vec q_vec fname_nsr

%--- Supress the full rank warning with cancorr
warning('off', 'stats:canoncorr:NotFullRank');

%--- For synthetic Dataset: Generate latent variables for predictors
% %- Number of samples (m), predictors (p) and responses (q)
% m      = 40;
% p      = 20;
% q      = 20;
% %- Noise level
% w_amp  = 0.5;
% %- Define LVs in each space
% d_XY   = 2;                % co-varying dimensions shared by X and Y
% d_Y    = d_XY*8;           % co-varying dimensions in Y only
% d_X    = d_XY*8;           % co-varying dimensions in X only
% d      = d_X + d_XY + d_Y; % total dimensionality of fluctuations
% Dx     = d_X(iex)+d_XY;
% idx_X  = (1:d_X);
% idx_XY = d_X + (1:d_XY);
% idx_Y  = d_X + d_XY + (1:d_Y);

%--- For real-world examples
% %- Corn dataset
% load('data/dataset_corn.mat')
% Dx = 4;         % Maxium LVs in X for corn dataset
%- pulp dataset
% load('data/milldata_kvarnsveden64.mat')
% %idxX = [1 2 3 4 5 6 7 9 10 12 13 15 16];      % 8 11 14 and 17 excluded
% %idxX = [1 5 6 7 8 9 10 11 12 13 14 15 16 17]; % 2 3 and excluded
% %idxY = [1 2 3 4 5 6 7 8];                     % 9 excluded
% %X    = X2z(:,idxX);
% %Y    = X3z(:,idxY);
% X     = X2z;
% Y     = X3z;
% Dx    = 5;             % [min=5, max=9=q] LVs in X for pulp dataset
%- Gene dataset
load('data/genes.mat')
Dx = 10; %[min=2, max=7]

%- Number of samples (m), predictors (p) and responses (q)
[m, p] = size(X);
q      = size(Y, 2);

%- Set k and r
max_components = min(p, q);
%max_components= min(d_X(iex)+d_XY, d_Y(iey)+d_XY);
r              = max_components;
k              = Dx;
if (k>r); error('Error: Violated k <= min(p,q).');  end

% Reproducibility parameters: Repitions, MCreps, CV
n_random_seeds = 20;  % Random numbers
n_repetitions  = 50;  % Number of Monte Carlo repetitions over CV
CV             = 10;  % 10-folsd cross-validation

% Preallocate arrays to store MSEs and other metrics
MSEmox_all      = zeros(n_random_seeds, max_components);
MSEmox_kisl_all = zeros(n_random_seeds, max_components);
MSEpls_all      = zeros(n_random_seeds, max_components);
MSEcca_all      = zeros(n_random_seeds, max_components);
MSEols_all      = zeros(n_random_seeds, 1);

for seed = 1:n_random_seeds
    rng(seed); % Set seed for reproducibility

    % --- MOX Regression ---
    MSEmox      = zeros(max_components, 1);
    MSEmox_kisl = zeros(max_components, 1);
    for h = 1:max_components
        % if max_component = min(p,q), and min(p,q) > max(Dx,Dy)
        k = max(h, Dx);
        [~, ~, ~, ~, ~, ~, ~, MSEcv, Fmaxcv, ~, ~, ~] = moxregress(X, Y, ...
            k, h, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox(h) = MSEcv/q;
        [~, ~, ~, ~, ~, ~, ~, MSEcv, Fmaxcv, ~, ~, ~] = moxregress(X, Y, ...
            h, h, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox_kisl(h) = MSEcv/q;
    end
    % Store results
    MSEmox_all(seed, :)      = MSEmox';
    MSEmox_kisl_all(seed, :) = MSEmox_kisl';

    % --- PLS Regression ---
    MSEpls = zeros(max_components, 1);
    for l = 1:max_components
        % Use the built-in cross-validation function of plsregress
        [~, ~, ~, ~, ~, ~, MSEcv] = plsregress(X, Y, l, 'CV', CV, 'MCReps', n_repetitions);
        MSEpls(l) = MSEcv(2, end)/q;
    end
    % Store results
    MSEpls_all(seed, :) = MSEpls';

    % --- OLS Regression (No components to tune) ---
    MSEols = 0;
    cvp = cvpartition(m, 'KFold', CV);
    for i = 1:CV
        trainIdx = training(cvp, i);
        testIdx  = test(cvp, i);

        % OLS regression using training data
        X_train = X(trainIdx, :);
        Y_train = Y(trainIdx, :);
        B_ols   = X_train \ Y_train;  % OLS solution

        % Predict on the test set
        X_test     = X(testIdx, :);
        Y_test     = Y(testIdx, :);
        Y_pred_ols = X_test * B_ols;

        % Compute MSE for the test set
        residuals = Y_test - Y_pred_ols;
        mse_fold  = mean(residuals(:).^2);
        MSEols    = MSEols + mse_fold / CV;  % Average over the folds
    end
    % Store result
    MSEols_all(seed) = MSEols;

    % --- CCA Regression ---
    MSEcca = zeros(max_components, 1);
    for l = 1:max_components
        MSEcv_total = 0;
        for i = 1:CV
            % Get training and test indices
            trainIdx = training(cvp, i);
            testIdx  = test(cvp, i);

            % Split the data
            X_train = X(trainIdx, :);
            Y_train = Y(trainIdx, :);
            X_test  = X(testIdx, :);
            Y_test  = Y(testIdx, :);

            % Perform CCA on training data
            [A, B, ~] = canoncorr(X_train, Y_train);

            % Keep only the first 'l' components
            A = A(:, 1:l);
            B = B(:, 1:l);

            % Project training data onto canonical variates
            U_train = X_train * A;
            V_train = Y_train * B;

            % Fit regression model from U_train to V_train
            beta = (U_train \ V_train)';  % Transpose to get coefficients as columns

            % Predict canonical variates of Y from test X
            U_test = X_test * A;
            V_pred = (beta * U_test')';  % Predicted canonical variates of Y

            % Reconstruct Y from predicted canonical variates
            Y_pred = V_pred * B(:, 1:l)';  % V_pred * B' gives Y in original space

            % Compute MSE for the test set
            residuals = Y_test - Y_pred;
            mse_fold  = mean(residuals(:).^2);

            MSEcv_total = MSEcv_total + mse_fold / CV;
        end
        MSEcca(l) = MSEcv_total;
    end
    % Store results
    MSEcca_all(seed, :) = MSEcca';
    
    % print repetition seed
    if mod(seed,5) == 0
        fprintf('\n i_repetition_seed:  %4d', seed);
    end
end
fprintf('\n');

% Compute mean MSEs across all random seeds
MSEmox_mean	     = mean(MSEmox_all, 1);
MSEmox_kisl_mean = mean(MSEmox_kisl_all, 1);
MSEpls_mean      = mean(MSEpls_all, 1);
MSEcca_mean      = mean(MSEcca_all, 1);
MSEols_mean      = mean(MSEols_all);

% Compute standard deviations if desired
MSEmox_std = std(MSEmox_all, 0, 1);
MSEmox_kisl_std = std(MSEmox_kisl_all, 0, 1);
MSEpls_std = std(MSEpls_all, 0, 1);
MSEcca_std = std(MSEcca_all, 0, 1);
MSEols_std = std(MSEols_all);

% Prepare data for plotting OLS with error bars
% Replicate the mean and standard deviation across the x-axis
ols_x = [0, max_components]; % x-axis range
ols_y = [MSEols_mean, MSEols_mean]; % Mean MSE replicated
ols_y_upper = ols_y + MSEols_std; % Upper bound (mean + std)
ols_y_lower = ols_y - MSEols_std; % Lower bound (mean - std)

% Plot the average cross-validated MSE as a function of the number of components
figure(1); hold on;
% OLS (horizontal line)
errorbar(max_components*[-1 0.5 2], MSEols_mean*[1 1 1], MSEols_std*[1 1 1], ...
    '--c', 'LineWidth', 3,'DisplayName', 'OLS');
% CCA
errorbar(1:max_components, MSEcca_mean, MSEcca_std, '-', ...
    'color', [0.72,0.27,1], ...
    'LineWidth', 3, 'Marker', 's', 'MarkerSize', 14,  ...
    'MarkerFaceColor', 'w', 'DisplayName', 'CCA', 'LineWidth',3);
% PLS
errorbar(1:max_components, MSEpls_mean, MSEpls_std, '-r', ...
    'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'PLS');
% MOX
errorbar(1:max_components, MSEmox_kisl_mean, MSEmox_kisl_std, '-k', ...
    'LineWidth', 3, 'Marker', '^','MarkerSize', 14, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'MOX_{hh}');
errorbar(1:max_components, MSEmox_mean, MSEmox_std, '-b', ...
    'LineWidth', 3, 'Marker', 's', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'MOX_{kh}');
xlabel('Number of Components (h)');
ylabel('Cross-validated MSE');
legend('Location', 'northeast');
xlim([0 max_components]);
ylim([0.0 * min([MSEmox_mean, MSEpls_mean, MSEols_mean]), ...
    2.0 * max([MSEmox_mean, MSEpls_mean, MSEols_mean])]);
grid on; set(gca,'FontSize', 24);  set(gcf, 'Color', 'w');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
hold off;
saveas(gcf,[ '03_gene_MSE']); close all;
save('03_results_geneDataset_vishal.mat');
