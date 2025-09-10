function vif_values = computeVIF(X)
p = size(X,2);
vif_values = zeros(1, p);
for i = 1:p
    X_i = X(:, i);
    X_other = X(:, setdiff(1:p, i));
    b = regress(X_i, X_other);
    Xi_hat = X_other * b;
    R2 = 1 - sum((X_i - Xi_hat).^2) / sum((X_i - mean(X_i)).^2);
    vif_values(i) = 1 / (1 - R2);
end
end