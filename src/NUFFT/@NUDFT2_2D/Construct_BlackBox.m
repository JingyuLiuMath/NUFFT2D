function Construct_BlackBox(A, xy, min_points, target_rank, tol)
% Construct_BlackBox

% Jingyu Liu, January 27, 2025.

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    min_points (1, 1) double;
    target_rank (1, 1) double;
    tol (1, 1) double;
end

A.M_ = size(xy, 1);
nx = A.nx_;
ny = A.ny_;
M = A.M_;
N = A.N_;

A.AFinv_HSS_ = NUDFT2_2D_HSS(nx, ny);
eta = ceil(M / N);
[A.xy_perm_, A.omega_perm_] = A.AFinv_HSS_.BuildTree(xy, min_points, eta, 0);
[~, A.xy_inv_perm_] = sort(A.xy_perm_, "ascend");
[~, A.omega_inv_perm_] = sort(A.omega_perm_, "ascend");

[row_leaf_size, col_leaf_size] = A.AFinv_HSS_.MaxLeafSize();

xy_p = xy(A.xy_perm_, :);
op_tildeA = @(v) op_tildeA_fun(v, xy_p, nx, ny, A.omega_inv_perm_);
op_tildeA_star = @(v) op_tildeA_fun_star(v, xy_p, nx, ny, A.omega_perm_);
A.AFinv_HSS_.BlackBoxConstruct(op_tildeA, op_tildeA_star, ...
    row_leaf_size, col_leaf_size, target_rank, tol);

end

function f = op_tildeA_fun(f, xy, nx, ny, omega_inv_perm)

arguments (Input)
    f (:, :) double;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    omega_inv_perm (:, 1) double;
end

arguments (Output)
    f (:, :) double;
end

N = nx * ny;
f = f(omega_inv_perm, :);
f = reshape(f, nx, ny, []);
f = ifft2(f);
f = reshape(f, N, []);
f = MY_NUFFT2_2D(f, xy, nx, ny);

end

function f = op_tildeA_fun_star(f, xy, nx, ny, omega_perm)

arguments (Input)
    f (:, :) double;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    omega_perm (:, 1) double;
end

arguments (Output)
    f (:, :) double;
end

N = nx * ny;
omega_x = nx * xy(:, 1);
omega_y = ny * xy(:, 2);
omega = [omega_x, omega_y];
f = MY_NUFFT1_2D(f, -omega, nx, ny);
f = reshape(f, nx, ny, []);
f = fft2(f) / N;
f = reshape(f, N, []);
f = f(omega_perm, :);

end