function A = NUDFT_Matrix(x, omega)
% NUDFT_Matrix


arguments (Input)
    x (:, 1) double;
    omega (:, 1) double;
end


arguments (Output)
    A (:, :) double;
end

A = exp(-2 * pi * 1i * x * omega.');

end