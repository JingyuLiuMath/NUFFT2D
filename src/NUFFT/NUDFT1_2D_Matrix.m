function A = NUDFT1_2D_Matrix(omega, mx, my)
% NUDFT_Matrix


arguments (Input)
    omega (:, 2) double;
    mx (1, 1) double;
    my (1, 1) double;
end


arguments (Output)
    A (:, :) double;
end

N = size(omega, 1);
M = mx * my;

Ax = NUDFT1_Matrix(omega(:, 1), mx);
Ay = NUDFT1_Matrix(omega(:, 2), my);

A = zeros(M, N);
for cit = 1 : N
    A(:, cit) = kron(Ay(:, cit), Ax(:, cit));
end

end