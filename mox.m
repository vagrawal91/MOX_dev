function [P, D, Q, muX, muY, E, Fmax, A, B, W] = mox(X, Y, k, l)
% MOX Perform Multivariate Orthogonal Cross-covariance (MOX) Regression.
%
%   This function performs MOX regression, which identifies latent variables 
%   through the maximization of cross-covariance between predictor and response 
%   matrices while ensuring orthogonality of the projections. It is based on 
%   Singular Value Decomposition (SVD) of the cross-covariance matrix.
%
%   [P, D, Q, muX, muY, E, Fmax, A, B, W] = mox(X, Y, k, l) returns:
%   
%   Inputs:
%       X  - Predictor matrix (m-by-p), where m is the number of observations 
%            and p is the number of predictor variables.
%       Y  - Response matrix (m-by-q), where m is the number of observations 
%            and q is the number of response variables.
%       k  - Number of predictor latent variables to retain (must satisfy k <= p).
%       l  - Number of response latent variables to retain (must satisfy l <= q).
%
%   Outputs:
%       P      - Final orthogonal projection matrix for predictors (p-by-l).
%       D      - Diagonal matrix of singular values (l-by-l), representing the 
%       Q      - Final orthogonal projection matrix for responses (q-by-l).
%                strength of each latent variable.
%       muX    - Mean of the predictors (1-by-p vector).
%       muY    - Mean of the responses (1-by-q vector).
%       E      - Transformed residuals (m-by-l).
%       Fmax   - Sum of the largest l singular values, representing the maximum 
%                cross-covariance captured by the latent variables.
%       A      - Intermediate orthogonal projection matrix for predictors (p-by-k).
%       B      - Intermediate orthogonal projection matrix for responses (p-by-l).
%       W      - Regression matrix (k-by-l).
%
%   Notes:
%       - The matrices X and Y must have the same number of rows (observations).
%       - The function centers both X and Y before computing the cross-covariance matrix.
%       - The function uses Singular Value Decomposition (SVD) to identify latent variables.
%       - Ensure that k <= p and l <= min(k, q) for valid input.
%
%   See also: svd, pinv

  % Dimensions of the input matrices
  m = size(X, 1); % number of observations (rows)
  p = size(X, 2); % number of predictor variables (columns)
  q = size(Y, 2); % number of response variables (columns)

  % Check if the number of latent variables l is valid
  if (k > min(p,q))
    error('Violated k <= min(p,q).');
  end
  if (l > k)
    error('Violated l <= k.');
  end

  % Center the predictor matrix X
  muX = mean(X);           % Mean of predictors
  X = X - ones(m, 1)*muX;  % Center predictors

  % Center the response matrix Y
  muY = mean(Y);           % Mean of responses
  Y = Y - ones(m, 1)*muY;  % Center responses

  % Cross-covariance matrix between centered X and Y
  SigmaXY = X' * Y;

  % Singular Value Decomposition of the cross-covariance matrix
  [U, S, V] = svd(SigmaXY, 'econ');
  A = U(:, 1:k);  % First k left singular vectors (predictors)
  B = V(:, 1:l);  % First l right singular vectors (responses)
  s = diag(S);    % Singular values from SVD
  Fmax = sum(s(1:l));  % Maximum value of the objective function

  % Latent variables for predictors and responses
  T = X * A;  % Latent variables for X (T) %v [mxp] [pxk]
  R = Y * B;  % Latent variables for Y (R) %v 

  % Regression matrix for latent variables and residual calculation
  W = pinv(T' * T) * T' * R;  % Solve T*W ≈ R using pseudoinverse
  e = R - T * W;  % Residuals (Y approximation error)

  % Singular Value Decomposition of the regression matrix W
  [Mp, Dp, N] = svd(W);
  M = Mp(:,1:l);
  D = Dp(1:l,:);
  
  % Adjust orthogonal transformations for final projections
  P = A * M;  % Final projection matrix for predictors    %v [pxh] <- [pxk][kxh], Actually M is [kxk] but M = M(:,1:h), as more h rows are zeros in D
  Q = B * N;  % Final projection matrix for responses     %v [qxh]
  E = e * N;  % Transformed residuals                     %v [mxh]

end

