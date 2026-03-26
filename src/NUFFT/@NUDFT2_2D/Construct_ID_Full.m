function Construct_ID_Full(A, xy, min_points, rank_or_tol)

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

x_col_pos = (0 : (nx - 1))';
y_col_pos = (0 : (ny - 1))';
col_pos = TensorProduct2D(x_col_pos, y_col_pos);

xy_p = xy(A.xy_perm_, :);
col_pos_q = col_pos(A.omega_perm_, :);
A.AFinv_HSS_.Construct_ID_Full(xy_p, col_pos_q, rank_or_tol);

end