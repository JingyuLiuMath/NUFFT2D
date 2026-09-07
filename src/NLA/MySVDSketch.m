function [U, S, V] = MySVDSketch(A, rank_or_tol)
% MySVDSketch

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    U (:, :) double;
    S (:, :) double;
    V (:, :) double;
end

[m, n] = size(A);
if isempty(A)
    U = zeros(m, 0);
    S = zeros(0);
    V = zeros(n, 0);
    return;
end

[U, S, V] = svd(A, "econ");
sigma = diag(S);

if rank_or_tol >= 1
    threshold = max(m, n) * eps(sigma(1));
else
    threshold = rank_or_tol * sigma(1);
end

k = find(sigma <= threshold, 1) - 1;
if isempty(k)
    k = length(sigma);
end

if rank_or_tol >= 1
    k = min(k, rank_or_tol);
end

U = U(:, 1 : k);
S = S(1 : k, 1 : k);
V = V(:, 1 : k);

end
