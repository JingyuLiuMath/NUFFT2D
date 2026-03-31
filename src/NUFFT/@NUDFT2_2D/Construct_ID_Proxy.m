function Construct_ID_Proxy(A, xy, min_points, rank_or_tol)

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

A.AFinv_HSS_.Construct_ID_Proxy(rank_or_tol);

end