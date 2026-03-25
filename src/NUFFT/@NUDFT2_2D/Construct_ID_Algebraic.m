function Construct_ID_Algebraic(A, xy, min_points, rank_or_tol)
% Construct

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    min_points (1, 1) double;
    rank_or_tol (1, 1) double;
end

A.M_ = size(xy, 1);
nx = A.nx_;
ny = A.ny_;
M = A.M_;
N = A.N_;

A.AFinv_HSS_ = NUDFT2_2D_HSS(nx, ny);
eta = ceil(M / N);
[A.xy_perm_, A.omega_perm_] = A.AFinv_HSS_.BuildTree(xy, min_points, eta, 1);
[~, A.xy_inv_perm_] = sort(A.xy_perm_, "ascend");
[~, A.omega_inv_perm_] = sort(A.omega_perm_, "ascend");
A.AFinv_HSS_.AssignHSSRowColumn(A.xy_perm_, A.omega_perm_);

H = NUDFT2_2D_Matrix(xy, nx, ny);

% Fx = NUDFT0_Matrix(nx, nx);
% Fy = NUDFT0_Matrix(ny, ny);
% F = kron(Fy, Fx);
% H = H / F;

H = reshape(H, M, nx, ny);
H = ifft(H, nx, 2);
H = ifft(H, ny, 3);
H = reshape(H, M, N);

A.AFinv_HSS_.Construct_ID_Algebraic(H, rank_or_tol);

end