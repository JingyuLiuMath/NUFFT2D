function [sk, re, T, k] = ID(A, rank_or_tol)
% ID Interpolative decomposition. A(:, re) = A(:, sk) * T

% Jingyu Liu, December 11, 2023.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    sk (1, :) double;
    re (1, :) double;
    T (:, :) double;
    k (1, 1) double;
end

[~, R, p, k] = MyQRSketch(A, rank_or_tol);

sk = p(1 : k);
re = p((k + 1) : end);
T = R(1 : k, 1 : k) \ R(1 : k, (k + 1) : size(R, 2));

% Exact ID.
% [~, n] = size(A);
% sk = 1 : n;
% re = [];
% T = zeros(n, 0);
% k = n;

% rel_err = norm(A(:, re) - A(:, sk) * T, "fro") / norm(A, "fro");
% fprintf("rel_err: %.1e\n", rel_err);

end