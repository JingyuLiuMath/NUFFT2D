function [Q, R, p, k] = MyQRSketch(A, rank_or_tol)
% MyQRSketch For a m-by-n matrix A, compute A(:, p) = Q * R where Q is
% m-by-k and R is k-by-n.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    Q (:, :) double;
    R (:, :) double;
    p (1, :) double;
    k (1, 1) double;
end

[m, n] = size(A);
if isempty(A)
    Q = zeros(m, 0);
    R = zeros(0, n);
    p = 1 : n;
    k = 0;
    return;
end

[Q, R, p] = qr(A, "econ", "vector");
if isvector(R)
    d = abs(R(1));
else
    d = abs(diag(R));
end

if rank_or_tol >= 1
    % The requested rank is an upper bound on the numerical rank.
    threshold = max(m, n) * eps(d(1));
else
    threshold = 0.1 * rank_or_tol * d(1);
end

k = find(d <= threshold, 1) - 1;
if isempty(k)
    k = length(d);
end

if rank_or_tol >= 1
    k = min(k, rank_or_tol);
end

Q = Q(:, 1 : k);
R = R(1 : k, :);

if nargout < 3
    [~, p_inv] = sort(p);
    R = R(:, p_inv);
end

end