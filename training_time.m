% Assuming q < m
% For OLS, dominated by
%   QR decomposition take O(mp²)
%   Solve triangular system O(p²q)
% For CCA, dominated by
%   Cross-covariance matrix O(mpq)
%   SVD take O(p²q) + O(pq²)
% For PLS with l latent variables, the cost is
%   Cross-covariance matrix O(mpq)
%   SVD take O(lp²q) + O(lpq²)
% For MOX, assuming l < r, dominated by
%   Cross-covariance matrix O(mpq)
%   SVD take O(p²q) + O(pq²)


% Number of random seeds (repetitions)
n_random_seeds = 20;

% Number of samples
m = 100;

% Number of predictor variables (p) and response variables (q)
p = 100;
q = 100;

% Noise level
w_amp = 0.5;

% Set the maximum number of components based on the formula
max_components = 6; % min(p, q);

% Preallocate arrays to store computational times
time_mox = zeros(n_random_seeds, 1);
time_pls = zeros(n_random_seeds, 1);
time_ols = zeros(n_random_seeds, 1);
time_cca = zeros(n_random_seeds, 1);

for seed = 1:n_random_seeds
    rng(seed); % Set seed for reproducibility

    % Generate latent variables for predictors
    d_X  = 2;   % co-varying dimensions in X only
    d_XY = 3;   % co-varying dimensions shared by X and Y
    d_Y  = 2;   % co-varying dimensions in Y only
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

    % --- Measure computational time for MOX Regression ---
    tic;
    for l = 1:max_components
        % Use MOX regression function
	[~, ~, ~] = mox(X, Y, l);
    end
    time_mox(seed) = toc;

    % --- Measure computational time for PLS Regression ---
    tic;
    for l = 1:max_components
        % Use PLS regression
        [~, ~] = plsregress(X, Y, l);
    end
    time_pls(seed) = toc;

    % --- Measure computational time for OLS Regression ---
    tic;
    B_ols = X \ Y;  % OLS solution
    time_ols(seed) = toc;

    % --- Measure computational time for CCA Regression ---
    tic;
    [A, B] = canoncorr(X, Y);  % Perform CCA
    for l = 1:max_components
        A_l = A(:, 1:l);
        B_l = B(:, 1:l);
    end
    time_cca(seed) = toc;
end

% Compute mean computational times across all random seeds
mean_time_mox = mean(time_mox);
mean_time_pls = mean(time_pls);
mean_time_ols = mean(time_ols);
mean_time_cca = mean(time_cca);

% Display the results
fprintf('Mean Computational Time (MOX): %.7f seconds\n', mean_time_mox);
fprintf('Mean Computational Time (PLS): %.7f seconds\n', mean_time_pls);
fprintf('Mean Computational Time (OLS): %.7f seconds\n', mean_time_ols);
fprintf('Mean Computational Time (CCA): %.7f seconds\n', mean_time_cca);
