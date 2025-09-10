function [T_p, P_p, W_p, C_p, T_o, P_o, W_o, B_p, X_residual, Y_residual] = opls(X, Y, num_predictive, num_orthogonal)
% OPLS - Orthogonal Partial Least Squares Regression
%
% Syntax:
%   [T_p, P_p, W_p, C_p, T_o, P_o, W_o, B_p, X_residual, Y_residual] = opls(X, Y, num_predictive, num_orthogonal)
%
% Inputs:
%   X               - Predictor matrix (m x p)
%   Y               - Response matrix (m x q)
%   num_predictive  - Number of predictive components to extract
%   num_orthogonal  - Number of orthogonal components to extract
%
% Outputs:
%   T_p        - Predictive scores (m x num_predictive)
%   P_p        - Predictive loadings (p x num_predictive)
%   W_p        - Predictive weights (p x num_predictive)
%   C_p        - Predictive Y-loadings (q x num_predictive)
%   T_o        - Orthogonal scores (m x num_orthogonal)
%   P_o        - Orthogonal loadings (p x num_orthogonal)
%   W_o        - Orthogonal weights (p x num_orthogonal)
%   B_p        - Regression coefficients for predictive components (p x q)
%   X_residual - Residual matrix of X after removing components
%   Y_residual - Residual matrix of Y after removing predictive components
%
% Example:
%   [T_p, P_p, W_p, C_p, T_o, P_o, W_o, B_p, X_residual, Y_residual] = opls(X, Y, 1, 1);

% Ensure that X and Y are centered
X = X - mean(X);
Y = Y - mean(Y);

[m, p] = size(X);
[~, q] = size(Y);

% Initialize residual matrices
E = X;
F = Y;

% Initialize matrices to store components
T_p = [];
P_p = [];
W_p = [];
C_p = [];
T_o = [];
P_o = [];
W_o = [];

% Extract predictive components
for i = 1:num_predictive
    % Calculate weights w_p
    w_p = E' * F / (F' * F);
    w_p = w_p / norm(w_p);  % Normalize w_p
    W_p = [W_p, w_p];
    
    % Calculate predictive score t_p
    t_p = E * w_p;
    T_p = [T_p, t_p];
    
    % Calculate predictive loading p_p
    p_p = E' * t_p / (t_p' * t_p);
    P_p = [P_p, p_p];
    
    % Calculate Y-loading c_p
    c_p = F' * t_p / (t_p' * t_p);
    C_p = [C_p, c_p];
    
    % Deflate E and F by removing predictive component
    E = E - t_p * p_p';
    F = F - t_p * c_p';
    
    % Extract orthogonal components
    for j = 1:num_orthogonal
        % Orthogonal weights w_o
        w_o = E' * (E * w_p);
        w_o = w_o - (w_p' * w_o) * w_p;  % Make w_o orthogonal to w_p
        w_o = w_o / norm(w_o);  % Normalize w_o
        W_o = [W_o, w_o];
        
        % Orthogonal score t_o
        t_o = E * w_o;
        T_o = [T_o, t_o];
        
        % Orthogonal loading p_o
        p_o = E' * t_o / (t_o' * t_o);
        P_o = [P_o, p_o];
        
        % Deflate E by removing orthogonal component
        E = E - t_o * p_o';
    end
end

% Regression coefficients for predictive components
B_p = W_p * inv(P_p' * W_p) * C_p';

% Residual matrices
X_residual = E;
Y_residual = F;

end
