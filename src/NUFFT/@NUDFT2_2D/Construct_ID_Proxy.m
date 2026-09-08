function Construct_ID_Proxy(A, xy, N_leaf, rank_or_tol)

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    N_leaf (1, 1) double;
    rank_or_tol (1, 1) double;
end

A.Factor_ = struct();
A.M_ = size(xy, 1);
nx = A.nx_;
ny = A.ny_;

A.AFinv_HSS_ = NUDFT2_2D_HSS(nx, ny);
[A.xy_perm_, A.omega_perm_] = A.AFinv_HSS_.Construct_ID_Proxy( ...
    xy, nx, ny, N_leaf, rank_or_tol);
[~, A.xy_inv_perm_] = sort(A.xy_perm_, "ascend");
[~, A.omega_inv_perm_] = sort(A.omega_perm_, "ascend");

end
