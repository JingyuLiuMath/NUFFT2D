function tildeA = NUFFT2_Kernel_Real(x, y, N)
% NUFFT2_Kernel_Real

arguments (Input)
    x (:, 1) double;
    y (:, 1) double;
    N (1, 1) double;
end

arguments (Output)
    tildeA (:, :) double;
end

phi = -2 * pi * 1i * (x - y.');
tildeA = expm1(-2 * pi * 1i * x * N) ./ expm1(phi) / N;

end