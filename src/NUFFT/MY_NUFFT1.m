function f = MY_NUFFT1(u, omega, M, method, acc)
% MY_NUFFT1

% f[k] = \sum_{j = 0}^{M - 1} exp(-2 * pi * 1i * (j / M) * omega[k]) * u[j]

arguments (Input)
    u (:, :) double;
    omega (:, 1) double;
    M (1, 1) double;
    method string = "finufft";
    acc (1, 1) double = 1e-15;
end


arguments (Output)
    f (:, :) double;
end

switch method
    case "direct"
        f = NUDFT1_Matrix(omega, M) * u;
    case "matlab"
        x = (0 : (M - 1))' / M;
        f = nufft(u, omega, x);
    case "finufft"
        f = finufft1d1(2 * pi * omega / M, ...
            exp(-2 * pi * 1i * floor(M / 2) * omega / M) .* u, ...
            -1, acc, M);
end

end