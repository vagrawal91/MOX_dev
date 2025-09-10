clear all;
close all;

load data/genes.mat

tic

% Centering and standardization
X = X - mean(X);
Y = zscore(Y);

% Sortera ut de största signalerna
rms = sqrt(sum(X.^2, 2)/size(X, 2)) * ones(1, size(X, 2));
bigs = any(abs(X) > 10*rms, 1);
X = X(:,bigs);


% Number of samples (m), predictor variables (p) and response variables (q)
[m, p] = size(X);
q = size(Y, 2);

% Set the maximum number of components based on the formula
r = min(p, q);

% Preallocate arrays to store MSEs and other metrics
n_random_seeds = 20;
MSEmox_rl_all = zeros(n_random_seeds, r);
MSEpls_all = zeros(n_random_seeds, r);

% Cross-validation settings
CV = 10; % 10-fold cross-validation
n_repetitions = 20; % Number of Monte Carlo repetitions within CV

for seed = 1:n_random_seeds
    seed
    rng(seed); % Set seed for reproducibility

    % --- MOX Regression ---
    MSEmox = zeros(r, 1);
    for l = 1:r
	[~, ~, ~, ~, ~, ~, ~, MSEcv, Fmaxcv, ~, ~, ~] = moxregress(X, Y, r, l, 'CV', CV, 'MCReps', n_repetitions);
        MSEmox(l) = MSEcv / q;
    end
    % Store results
    MSEmox_rl_all(seed, :) = MSEmox';

    % --- PLS Regression ---
    MSEpls = zeros(r, 1);
    for l = 1:r
        % Use the built-in cross-validation function of plsregress
        [~, ~, ~, ~, ~, ~, MSEcv] = plsregress(X, Y, l, 'CV', CV, 'MCReps', n_repetitions);
        MSEpls(l) = MSEcv(2, end) / q;
    end
    % Store results
    MSEpls_all(seed, :) = MSEpls';

end

% Compute mean MSEs across all random seeds
MSEmox_rl_mean = mean(MSEmox_rl_all, 1);
MSEpls_mean = mean(MSEpls_all, 1);

% Compute standard deviations if desired
MSEmox_rl_std = std(MSEmox_rl_all, 0, 1);
MSEpls_std = std(MSEpls_all, 0, 1);

save('result_gene_dataset.mat', 'MSEmox_rl_mean', 'MSEpls_mean', 'MSEmox_rl_std', 'MSEpls_std', 'p', 'q', 'm', 'r', '-nocompression');

figure(1)
plot(1:r, MSEpls_mean, 1:r, MSEmox_rl_mean)
legend('PLS','MOX_{\itrl}')		      

toc

