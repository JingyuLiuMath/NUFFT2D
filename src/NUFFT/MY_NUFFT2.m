function f = MY_NUFFT2(u, x, N, method, acc)
% MY_NUFFT2

% f[j] = \sum_{k = 0}^{N - 1} exp(-2 * pi * 1i * k * x[j]) u[k]
% for 0 <= j < M.
% N = size(u, 1).
% M = size(f, 1) = length(x).

% Jingyu Liu, November 17, 2024.

arguments (Input)
    u (:, :) double;
    x (:, 1) double;
    N (1, 1) double;
    method string = "finufft";
    acc (1, 1) double = 1e-15;
end


arguments (Output)
    f (:, :) double;
end

switch method
    case "direct"
        f = NUDFT2_Matrix(x, N) * u;
    case "matlab"
        f = nufft(u, [], x);
    case "finufft"
        f = exp(-2 * pi * 1i * floor(N / 2) * x) ...
            .* finufft1d2(2 * pi * x, -1, acc, u);
end

end