function Construct_ID_Proxy(A, xy, N_leaf, tol)

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    N_leaf (1, 1) double;
    tol (1, 1) double;
end

A.Factor_ = struct();
A.M_ = size(xy, 1);
nx = A.nx_;
ny = A.ny_;

A.AFinv_HSS_ = NUDFT2_2D_HSS(nx, ny);
[A.xy_perm_, A.omega_perm_] = A.AFinv_HSS_.BuildTree(xy, nx, ny, N_leaf);
A.AFinv_HSS_.Construct_ID_Proxy(xy(A.xy_perm_, :), nx, ny, tol);
A.xy_inv_perm_ = zeros(size(A.xy_perm_));
A.xy_inv_perm_(A.xy_perm_) = (1 : size(A.xy_perm_, 1))';
A.omega_inv_perm_ = zeros(size(A.omega_perm_));
A.omega_inv_perm_(A.omega_perm_) = (1 : size(A.omega_perm_, 1))';

end
