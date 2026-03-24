function A = NUDFT2_2D_Matrix(xy, nx, ny)
% NUDFT_Matrix


arguments (Input)
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
end


arguments (Output)
    A (:, :) double;
end

M = size(xy, 1);
N = nx * ny;

Ax = NUDFT2_Matrix(xy(:, 1), nx);
Ay = NUDFT2_Matrix(xy(:, 2), ny);

A = zeros(M, N);
for rit = 1 : M
    A(rit, :) = kron(Ay(rit, :), Ax(rit, :));
end

end