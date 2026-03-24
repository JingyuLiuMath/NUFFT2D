function [Q, k] = ColBasis(B, k, tol)
% ColBasis

% Jingyu Liu, December 4, 2024.

arguments (Input)
    B (:, :) double;
    k (1, 1) double;
    tol (1, 1) double;
end

arguments (Output)
    Q (:, :) double;
    k (1, 1) double;
end

[Q, R, ~] = qr(B, "econ", "vector");

k = min(find(abs(diag(R)) >= tol * 1e-1 * max(abs(R(1, 1)), 1), 1, "last"), k);

Q = Q(:, 1 : k); 

end