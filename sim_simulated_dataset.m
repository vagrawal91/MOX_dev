%clear all;
clear;%
close all;

% Number of random seeds (repetitions)
n_random_seeds = 20;

% Number of samples
m = 40;

% Number of predictor variables (p) and response variables (q)
p = 20;
q = 10;

% Noise level
w_amp = 0.5;

% Set the maximum number of components based on the formula
r = min(p, q);

% Preallocate arrays to store MSEs and other metrics
MSEmox_rl_all = zeros(n_random_seeds, r);
MSEmox_ll_all = zeros(n_random_seeds, r);
MSEpls_all = zeros(n_random_seeds, r);
MSEcca_all = zeros(n_random_seeds, r);
MSEols_all = zeros(n_random_seeds, 1);

% Cross-validation settings
CV = 10; % 10-fold cross-validation
n_repetitions = 50; % Number of Monte Carlo repetitions within CV

% Dimensionality of latent structures
d_X  = 2;   % co-varying dimensions in X only
d_XY = 3;   % co-varying dimensions shared by X and Y
d_Y  = 2;   % co-varying dimensions in Y only
d = d_X + d_XY + d_Y; % total dimensionality of fluctuations
idx_X  = (1:d_X);
idx_XY = d_X + (1:d_XY);
idx_Y  = d_X + d_XY + (1:d_Y);

for seed = 1:n_random_seeds
    rng(seed); % Set seed for reproducibility

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
    MSEmox = zeros(r, 1);
    MSEmox_ll = zeros(r, 1);
    for l = 1:r
	[~, ~, ~, ~, ~, ~, ~, MSEcv, ~, ~, ~, ~] = moxregress(X, Y, q, l, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox(l) = MSEcv / q;
	[~, ~, ~, ~, ~, ~, ~, MSEcv, ~, ~, ~, ~] = moxregress(X, Y, l, l, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox_ll(l) = MSEcv / q;
    end
    % Store results
    MSEmox_rl_all(seed, :) = MSEmox';
    MSEmox_ll_all(seed, :) = MSEmox_ll';

    % --- PLS Regression ---
    MSEpls = zeros(r, 1);
    for l = 1:r
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
    MSEcca = zeros(r, 1);
    for l = 1:r
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
MSEmox_rl_mean = mean(MSEmox_rl_all, 1);
MSEmox_ll_mean = mean(MSEmox_ll_all, 1);
MSEpls_mean = mean(MSEpls_all, 1);
MSEcca_mean = mean(MSEcca_all, 1);
MSEols_mean = mean(MSEols_all);

% Compute standard deviations if desired
MSEmox_rl_std = std(MSEmox_rl_all, 0, 1);
MSEmox_ll_std = std(MSEmox_ll_all, 0, 1);
MSEpls_std = std(MSEpls_all, 0, 1);
MSEcca_std = std(MSEcca_all, 0, 1);
MSEols_std = std(MSEols_all);

%save('result_simulated_dataset.mat', 'MSEmox_rl_mean', 'MSEmox_ll_mean', 'MSEpls_mean', 'MSEcca_mean', 'MSEols_mean', 'MSEmox_rl_std', 'MSEmox_ll_std', 'MSEpls_std', 'MSEcca_std', 'MSEols_std', 'p', 'q', 'm', 'r', '-nocompression');
