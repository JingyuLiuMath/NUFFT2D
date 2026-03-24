function A = NUDFT1_Matrix(omega, m)
% NUDFT_Matrix


arguments (Input)
    omega (:, 1) double;
    m (1, 1) double;
end


arguments (Output)
    A (:, :) double;
end

x = (0 : (m - 1))' / m;
A = exp(-2 * pi * 1i * x * omega.');

end