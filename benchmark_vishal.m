close all;
clearvars -except filename m p q w_amp n_random_seeds n_repetitions CV d_XY ...
    eps_x eps_y d_X d_XY d_Y d idx_X idx_XY idx_Y iex iey im iq ip ...
    m_vec p_vec q_vec

% Choose how X and Y are constructed. 1: in terms of Lambda, 2: Intutive
ini_opt        = 1;

% Set the maximum number of components based on the formula
max_components = min(p, q);
%max_components = d_Y(iey)+d_XY;
%max_components  = max(d_X(iex)+d_XY, d_Y(iey)+d_XY);

% Preallocate arrays to store MSEs and other metrics
MSEmox_all      = zeros(n_random_seeds, max_components);
MSEmox_kisl_all = zeros(n_random_seeds, max_components);
MSEpls_all      = zeros(n_random_seeds, max_components);
MSEcca_all      = zeros(n_random_seeds, max_components);
MSEols_all      = zeros(n_random_seeds, 1);

for seed = 1:n_random_seeds
    rng(seed); % Set seed for reproducibility
    
    if ini_opt == 1
        % Fluctuations of predictors
        t = randn(m, d); % each t(:, i) is a latent variable

        % Generate predictor matrix X using latent variables and random loadings
        xl = randn(p, d);  % latent
        xl(:, idx_Y) = 0;

        % Predictor matrix X
        X = zscore(t * xl' + randn(m, p) * w_amp);

        % Generate response variables
        yl = randn(q, d);
        yl(:, idx_X) = 0;

        % Response matrix Y
        Y = zscore(t * yl' + randn(m, q) * w_amp);

    elseif ini_opt == 2
        %-=-=-= new way of generating predictor and response arrays
        % Generate latent variables
        Z_X  = randn(m, d_X);       % latent for X only
        Z_Y  = randn(m, d_Y);       % latent for Y only
        Z_XY = randn(m, d_XY);      % shared latent
        
        % Generate random weights
        W_X    = 2 * rand(d_X, p);
        W_Y    = 2 * rand(d_Y, q);
        W_XY_X = 2 * rand(d_XY, p);
        W_XY_Y = 2 * rand(d_XY, q);
    
        % Generate Gaussian noise
        noise_X = w_amp * randn(m, p);
        noise_Y = w_amp * randn(m, q);
        
        % Construct observed data
        X = Z_X*W_X + Z_XY*W_XY_X + noise_X;
        Y = Z_Y*W_Y + Z_XY*W_XY_Y + noise_Y;
    
        % Standardization
        X = zscore(X);
        Y = zscore(Y);
    end
    
    % --- MOX Regression ---
    MSEmox      = zeros(max_components, 1);
    MSEmox_kisl = zeros(max_components, 1);
    for h = 1:max_components
        k               = max(h,d_X(iex)+d_XY);
        l               = max(h,d_Y(iey)+d_XY);
	    [~, ~, ~, ~, ~, ~, ~, MSEcv, Fmaxcv, ~, ~, ~] = moxregress(X, Y, ...
            k, l, h, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox(h) = MSEcv/q;
        [~, ~, ~, ~, ~, ~, ~, MSEcv, Fmaxcv, ~, ~, ~] = moxregress(X, Y, ...
            h, h, h, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox_kisl(h) = MSEcv / q;
    end
    % Store results
    MSEmox_all(seed, :)      = MSEmox';
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
		if mod(seed,5) == 0
			fprintf('i_repetition_seed:  %4d \n', seed);
		end
end

% Compute mean MSEs across all random seeds
MSEmox_mean			 = mean(MSEmox_all, 1);
MSEmox_kisl_mean = mean(MSEmox_kisl_all, 1);
MSEpls_mean      = mean(MSEpls_all, 1);
MSEcca_mean      = mean(MSEcca_all, 1);
MSEols_mean      = mean(MSEols_all);

% % Compute standard deviations if desired
% MSEmox_std = std(MSEmox_all, 0, 1);
% MSEmox_kisl_std = std(MSEmox_kisl_all, 0, 1);
% MSEpls_std = std(MSEpls_all, 0, 1);
% MSEcca_std = std(MSEcca_all, 0, 1);
% MSEols_std = std(MSEols_all);
% 
% % Prepare data for plotting OLS with error bars
% % Replicate the mean and standard deviation across the x-axis
% ols_x = [0, max_components]; % x-axis range
% ols_y = [MSEols_mean, MSEols_mean]; % Mean MSE replicated
% ols_y_upper = ols_y + MSEols_std; % Upper bound (mean + std)
% ols_y_lower = ols_y - MSEols_std; % Lower bound (mean - std)
% 
% % Plot the average cross-validated MSE as a function of the number of components
% figure(1);
% % OLS (horizontal line)
% errorbar(max_components*[-1 0.5 2], MSEols_mean*[1 1 1], MSEols_std*[1 1 1], '--', 'DisplayName', 'OLS', 'LineWidth',3);
% hold on;
% % CCA
% errorbar(1:max_components, MSEcca_mean, MSEcca_std, '-o', 'DisplayName', 'CCA', 'LineWidth',3);
% % PLS
% errorbar(1:max_components, MSEpls_mean, MSEpls_std, '-s', 'DisplayName', 'PLS', 'LineWidth',3);
% % MOX
% errorbar(1:max_components, MSEmox_kisl_mean, MSEmox_kisl_std, '-^', 'DisplayName', 'MOX_{ℓℓ}', 'LineWidth',3);
% errorbar(1:max_components, MSEmox_mean, MSEmox_std, '-d', 'DisplayName', 'MOX_{{\itr}ℓ}', 'LineWidth',3);
% xlabel('Number of Components (ℓ)');
% ylabel('Cross-validated MSE');
% legend('Location', 'northeast');
% xlim([0 max_components]);
% ylim([0.0 * min([MSEmox_mean, MSEpls_mean, MSEols_mean]), ...
%       2.0 * max([MSEmox_mean, MSEpls_mean, MSEols_mean])]);
% xticks([0:1:max_components]); yticks([0:0.1:2]);
% grid on; set(gca,'FontSize', 28);  set(gcf, 'Color', 'w');
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% hold off;

% MSE levels
%MSE_all = MSEvalues_all(d_X, d_XY, MSEmox_mean, MSEmox_kisl_mean, ...
%                    MSEpls_mean, MSEcca_mean, MSEols_mean);
%disp(MSE_all);
