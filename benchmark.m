clear;
close all;

%load dataset_tecator.mat
%load dataset_corn.mat

%load ../../modes/analysis/modes_kvarnsveden64.mat
%load('/Users/vishalagrawal/DriveD/mox_sim/sim_dimsame/data/milldata_kvarnsveden64.mat')

%idxX = [1 2 3 4 5 6 7 9 10 12 13 15 16]; % 8 11 14 and 17 excluded
%idxX = [1 5 6 7 8 9 10 11 12 13 14 15 16 17]; %v 2 3 and excluded
%idxY = [1 2 3 4 5 6 7 8]; % 9 excluded
%X = X2z(:,idxX);
%Y = X3z(:,idxY);

% Number of samples (m), predictor variables (p) and response variables (q)
%[m, p] = size(X);
%q      = size(Y, 2);

%v added
m = 40;
p = 20;
q = 40;

% Noise level
w_amp = 0.5;

% Set the maximum number of components based on the formula
max_components = min(p, q);

% Preallocate arrays to store MSEs and other metrics
n_random_seeds = 20;
MSEmox_all = zeros(n_random_seeds, max_components);
MSEmox_kisl_all = zeros(n_random_seeds, max_components);
MSEpls_all = zeros(n_random_seeds, max_components);
MSEcca_all = zeros(n_random_seeds, max_components);
MSEols_all = zeros(n_random_seeds, 1);

% Cross-validation settings
CV = 10; % 10-fold cross-validation
%n_repetitions = 20; % Number of Monte Carlo repetitions within CV
n_repetitions = 20; % Number of Monte Carlo repetitions over CV

for seed = 1:n_random_seeds
    rng(seed); % Set seed for reproducibility

    % Generate latent variables for predictors
    d_X  = 16;   % co-varying dimensions in X only
    d_XY = 2;   % co-varying dimensions shared by X and Y
    d_Y  = 16;   % co-varying dimensions in Y only
    d = d_X + d_XY + d_Y; % total dimensionality of fluctuations
    idx_X  = (1:d_X);
    idx_XY = d_X + (1:d_XY);
    idx_Y  = d_X + d_XY + (1:d_Y);

    % Fluctuations of predictors
    t = randn(m, d); % each t(:, i) is a latent variable

    % Generate predictor matrix X using latent variables and random loadings
    xl = randn(p, d);
    xl(:, idx_Y) = 0;

    % Predictor matrix X
    X = zscore(t * xl' + randn(m, p) * w_amp);

    % Generate response variables
    yl = randn(q, d);
    yl(:, idx_X) = 0;

    % Response matrix Y
    Y = zscore(t * yl' + randn(m, q) * w_amp);

    % --- MOX Regression ---
    MSEmox = zeros(max_components, 1);
    MSEmox_kisl = zeros(max_components, 1);
    for l = 1:max_components
	[~, ~, ~, ~, ~, ~, ~, MSEcv, Fmaxcv, ~, ~, ~] = moxregress(X, Y, q, l, 'CV', CV, 'MCReps', n_repetitions);
        %MSEmox(l) = MSEcv;
        MSEmox(l) = MSEcv/q;
	[~, ~, ~, ~, ~, ~, ~, MSEcv, Fmaxcv, ~, ~, ~] = moxregress(X, Y, l, l, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox_kisl(l) = MSEcv / q;
    end
    % Store results
    MSEmox_all(seed, :) = MSEmox';
    MSEmox_kisl_all(seed, :) = MSEmox_kisl';

    % --- PLS Regression ---
    MSEpls = zeros(max_components, 1);
    for l = 1:max_components
        % Use the built-in cross-validation function of plsregress
        [~, ~, ~, ~, ~, ~, MSEcv] = plsregress(X, Y, l, 'CV', CV, 'MCReps', n_repetitions);
        MSEpls(l) = MSEcv(2, end) / q;
    end
    % Store results
    MSEpls_all(seed, :) = MSEpls';

    % --- OLS Regression (No components to tune) ---
    MSEols = 0;
    cvp = cvpartition(m, 'KFold', CV);
    for i = 1:CV
        trainIdx = training(cvp, i);
        testIdx = test(cvp, i);

        % OLS regression using training data
        X_train = X(trainIdx, :);
        Y_train = Y(trainIdx, :);
        B_ols = X_train \ Y_train;  % OLS solution

        % Predict on the test set
        X_test = X(testIdx, :);
        Y_test = Y(testIdx, :);
        Y_pred_ols = X_test * B_ols;

        % Compute MSE for the test set
        residuals = Y_test - Y_pred_ols;
        mse_fold = mean(residuals(:).^2);

        MSEols = MSEols + mse_fold / CV;  % Average over the folds
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
            testIdx = test(cvp, i);

            % Split the data
            X_train = X(trainIdx, :);
            Y_train = Y(trainIdx, :);
            X_test = X(testIdx, :);
            Y_test = Y(testIdx, :);

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
            mse_fold = mean(residuals(:).^2);

            MSEcv_total = MSEcv_total + mse_fold / CV;
        end
        MSEcca(l) = MSEcv_total;
    end
    % Store results
    MSEcca_all(seed, :) = MSEcca';
end

% Compute mean MSEs across all random seeds
MSEmox_mean = mean(MSEmox_all, 1);
MSEmox_kisl_mean = mean(MSEmox_kisl_all, 1);
MSEpls_mean = mean(MSEpls_all, 1);
MSEcca_mean = mean(MSEcca_all, 1);
MSEols_mean = mean(MSEols_all);

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
figure(1);
% OLS (horizontal line)
errorbar(max_components*[-1 0.5 2], MSEols_mean*[1 1 1], MSEols_std*[1 1 1], '--', 'DisplayName', 'OLS');
hold on;
% CCA
errorbar(1:max_components, MSEcca_mean, MSEcca_std, '-o', 'DisplayName', 'CCA');
% PLS
errorbar(1:max_components, MSEpls_mean, MSEpls_std, '-s', 'DisplayName', 'PLS');
% MOX
errorbar(1:max_components, MSEmox_kisl_mean, MSEmox_kisl_std, '-^', 'DisplayName', 'MOX_{ℓℓ}');
errorbar(1:max_components, MSEmox_mean, MSEmox_std, '-d', 'DisplayName', 'MOX_{{\itr}ℓ}');
hold off;
xlabel('Number of Components (ℓ)');
ylabel('Cross-validated MSE');
legend('Location', 'northeast');
xlim([0 max_components]);
ylim([0.0 * min([MSEmox_mean, MSEpls_mean, MSEols_mean]), ...
      2.0 * max([MSEmox_mean, MSEpls_mean, MSEols_mean])]);