function [optimal_k, cumulative_variance] = determine_optimal_svd_components(X, varargin)
    % DETERMINE_OPTIMAL_SVD_COMPONENTS - Find optimal number of latent variables using SVD
    %
    % Determines the minimum number of SVD components needed to capture 98%
    % percentage of variance in a given predictor/response arrays.
    %
    % Inputs:
    %   X - Data matrix (samples x variables), must be numeric and 2D
    %
    % Optional Name-Value Pairs:
    %   'variance_threshold'  - Cumulative variance threshold 0-1 (default: 0.99)
    %   'max_components'      - Maximum components to consider (default: min(size(X))-1)
    %
    % Outputs:
    %   optimal_k - Minimum number of components to reach variance threshold
    
    % Parse and validate inputs
    p = inputParser;
    addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'variance_threshold', 0.98, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    addParameter(p, 'max_components', min(size(X))-1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, X, varargin{:});
    
    % Extract parameters
    var_threshold = p.Results.variance_threshold;
    max_k         = min(p.Results.max_components, min(size(X))-1);
    
    % Validate data
    if any(isnan(X(:))) || any(isinf(X(:)))
        error('X contains NaN or Inf values');
    end
    
    [n_samples, n_vars] = size(X);
    if n_samples < 2 || n_vars < 1
        error('X must have at least 2 samples and 1 variable');
    end
    
    % Center the data (critical for meaningful variance decomposition)
    X_centered = X - mean(X, 1);
    
    % Perform economical SVD
    [~, S, ~] = svd(X_centered, 'econ');
    singular_values = diag(S);
    
    % Limit to max_components
    n_available     = min(length(singular_values), max_k);
    singular_values = singular_values(1:n_available);
    
    % Calculate variance metrics
    variance_per_component = singular_values.^2;
    total_variance         = sum(variance_per_component);
    explained_variance     = variance_per_component / total_variance;
    cumulative_variance    = cumsum(explained_variance);
    
    % Find optimal number of components
    idx_threshold = find(cumulative_variance >= var_threshold, 1, 'first');
    
    if isempty(idx_threshold)
        % Threshold not reached - use all available components and warn
        warning('determine_optimal_svd_components: Threshold Not Met')
        idx_threshold = find(cumulative_variance >= 0.95, 1, 'first');
        optimal_k     = idx_threshold;
    else
        optimal_k = idx_threshold;
    end
    end