function A = NUDFT0_Matrix(m, n)
% NUDFT_Matrix


arguments (Input)
    m (1, 1) double;
    n (1, 1) double;
end


arguments (Output)
    A (:, :) double;
end

x = (0 : (m - 1))' / m;
omega = (0 : (n - 1))';
A = exp(-2 * pi * 1i * x * omega.');

end