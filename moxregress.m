function [P, D, Q, muX, muY, E, Fmax, MSEcv, Fmaxcv, A, B, W] = moxregress(X, Y, k, h, varargin)
% MOXREGRESS Perform MOX regression with optional cross-validation and multiple restarts.
%
%   [P, D, Q, muX, muY, E, Fmax, MSEcv, Fmaxcv, A, B, W] = moxregress(X, Y, k, l, 'CV', j, 'MCReps', n)
%   applies the MOX regression to the data in X and Y with l latent variables.
%
%   Inputs:
%       X  - Predictor matrix (m-by-p), where m is the number of observations 
%            and p is the number of predictor variables.
%       Y  - Response matrix (m-by-q), where m is the number of observations 
%            and q is the number of response variables.
%       k  - Number of predictor latent variables (must satisfy k <= p).
%       l  - Number of response latent variables to retain (must satisfy l <= q).
%
%   Optional parameters:
%     'CV'     - An integer number specifying j-fold cross-validation. Default is 1 (no CV).
%     'MCReps' - Number of Monte Carlo repetitions for cross-validation. Default is 1 (no restarts).
%
%   Outputs:
%     P       - Final orthogonal projection matrix for predictors (p-by-l).
%     D       - Diagonal matrix of singular values (l-by-l), representing the 
%     Q       - Final orthogonal projection matrix for responses (q-by-l).
%     muX     - Mean of the predictors (1-by-p vector).
%     muY     - Mean of the responses (1-by-q vector).
%     E       - Transformed residuals (m-by-l).
%     Fmax    - Maximum cross-covariance captured by latent variables.
%     MSEcv   - Cross-validated Mean Squared Error of residuals.
%     Fmaxcv  - Cross-validated Fmax value, computed based on the test set.
%     A       - Orthogonal projection matrix for predictors (p-by-k).
%     B       - Orthogonal projection matrix for responses (q-by-l).
%     W       - Regression matrix between latent variables (k-by-l).
%
%   See also: mox, plsregress, svd

% Default values for optional parameters
CV = 1;       % Default: no cross-validation
MCReps = 1;   % Default: no multiple restarts

% Parse optional parameters
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'CV')
        CV = varargin{i+1};  % Number of cross-validation folds
    elseif strcmpi(varargin{i}, 'MCReps')
        MCReps = varargin{i+1};  % Number of Monte Carlo repetitions
    end
end

% Initialize variables for cross-validation
[m, p] = size(X);
q = size(Y, 2);

MSEcv = zeros(MCReps, 1);
Fmaxcv = zeros(MCReps, 1);

% Create cross-validation partition if CV > 1
if CV > 1
    cvp = cvpartition(m, 'KFold', CV);
else
    cvp = cvpartition(m, 'Resubstitution');
end

% Cross-validation loop with MCReps repetitions
for rep = 1:MCReps
    mse_sum = 0;  % Sum of MSE  across all folds
    fmax_sum = 0; % Sum of Fmax across all folds

    % Loop over cross-validation folds
    for fold = 1:cvp.NumTestSets
        trainIdx = training(cvp, fold);  % Training indices
        testIdx  = test(cvp, fold);      % Test indices
        
        % Perform MOX on training data
        [P_train, D_train, Q_train, muX_train, muY_train, ~, ~, A_train, ...
            B_train, W_train] = mox(X(trainIdx, :), Y(trainIdx, :), k, l, h);
        
        % Center the test data using the means from the training set
        X_test = X(testIdx, :) - muX_train;
        Y_test = Y(testIdx, :) - muY_train;
        
        % Predict the responses for the test data
        Y_pred = X_test * P_train * D_train * Q_train';
        
        % Compute residuals for test data and MSE
        residuals = Y_test - Y_pred;
        mse_fold  = mean(residuals(:).^2);
        mse_sum   = mse_sum + mse_fold;
        
        % Calculate cross-covariance for test data
        Sigma_test = X_test' * Y_test;

        % Compute Fmax for the test set
        s_fold    = diag(A_train'*Sigma_test*B_train);  % Sum of the first l singular values for the test data
				Fmax_fold = sum(s_fold(1:min(k,h)));
        fmax_sum  = fmax_sum + Fmax_fold;
    end
    
    % Average MSE and Fmax over all folds
    MSEcv(rep) = q * mse_sum / CV; % This is to make MSE behave as in PLS regress, which sums across columns instead of taking average
    Fmaxcv(rep) = fmax_sum / CV;
end

% Compute final cross-validated metrics
MSEcv = mean(MSEcv);  % Average over MCReps
Fmaxcv = mean(Fmaxcv);  % Average over MCReps

% Perform final MOX on the entire dataset
[P, D, Q, muX, muY, E, Fmax, A, B, W] = mox(X, Y, k, h);
end
