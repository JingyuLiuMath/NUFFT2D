function f = MY_NUFFT(u, x, omega, method, acc)
% MY_NUFFT

% f[j] = \sum_{k = 0}^{N - 1} exp(-2 * pi * 1i * x[j] * omega[k]) u[k]
% for 0 <= j < M.
% N = size(u, 1) = length(omega).
% M = size(f, 1) = length(x).

% Jingyu Liu, November 17, 2024.

arguments (Input)
    u (:, :) double;
    x (:, 1) double = (0 : (size(u, 1) - 1))' / size(u, 1);
    omega (:, 1) double = (0 : (size(u, 1) - 1))';
    method string = "finufft";
    acc (1, 1) double = 1e-15;
end


arguments (Output)
    f (:, :) double;
end

switch method
    case "direct"
        f = NUDFT_Matrix(x, omega) * u;
    case "matlab"
        f = nufft(u, omega, x);
    case "finufft"
        f = finufft1d3(omega, u, -1, acc, 2 * pi * x);
end

end