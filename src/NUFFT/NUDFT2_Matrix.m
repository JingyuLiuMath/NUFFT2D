function A = NUDFT2_Matrix(x, n)
% NUDFT_Matrix


arguments (Input)
    x (:, 1) double;
    n (1, 1) double;
end


arguments (Output)
    A (:, :) double;
end

omega = (0 : (n - 1))';
A = exp(-2 * pi * 1i * x * omega.');

end