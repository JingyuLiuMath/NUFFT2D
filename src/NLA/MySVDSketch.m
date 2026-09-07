function [U, S, V] = MySVDSketch(A, rank_or_tol)
% MySVDSketch

% Jingyu Liu, May 9, 2024.

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

if rank_or_tol >= 1
    [U, S, V] = svds(A, rank_or_tol);
else
    if isempty(A)
        U = zeros(size(A, 1), 0);
        S = zeros(0);
        V = zeros(size(A, 2), 0);
        return;
    end
    [U, S, V] = svd(A, "econ");
    sigma = diag(S);
    k = find(sigma >= rank_or_tol * sigma(1) & sigma > 0, 1, "last");
    if isempty(k)
        k = 0;
    end
    U = U(:, 1 : k);
    S = S(1 : k, 1 : k);
    V = V(:, 1 : k);
end

end
