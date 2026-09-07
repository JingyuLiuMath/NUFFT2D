function [sk, re, T, k] = ID(A, rank_or_tol)
% ID Interpolative decomposition. A(:, re) = A(:, sk) * T

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
if k == 0
    T = zeros(0, numel(re));
else
    T = R(:, 1 : k) \ R(:, (k + 1) : end);
end

% rel_err = norm(A(:, re) - A(:, sk) * T, "fro") / norm(A, "fro");
% fprintf("rel_err: %.1e\n", rel_err);

end