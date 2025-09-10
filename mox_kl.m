function [P, D, Q] = mox_kl(X, Y, k, l)
  [U, ~, V] = svd(X' * Y, "econ");
  A = U(:, 1:k);
  B = V(:, 1:l);
  [M, D, N] = svd((X * A) \ (Y * B), "econ");
  P = A * M;
  Q = B * N;
end
