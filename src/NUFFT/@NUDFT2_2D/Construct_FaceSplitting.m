function Construct_FaceSplitting(A, xy, n_leaf, tol_hss)

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    n_leaf (1, 1) double;
    tol_hss (1, 1) double;
end

A.Factor_ = struct();
A.M_ = size(xy, 1);
nx = A.nx_;
ny = A.ny_;

Ax = NUDFT2_HSS(nx);
px = Ax.BuildTree(xy(:, 1), n_leaf, 2);
Ax.Construct_fADI(tol_hss / 10);
Ay = NUDFT2_HSS(ny);
py = Ay.BuildTree(xy(:, 2), n_leaf, 2);
Ay.Construct_fADI(tol_hss / 10);

A.AFinv_HSS_ = NUDFT2_2D_HSS(nx, ny);
[A.xy_perm_, A.omega_perm_] = ...
    A.AFinv_HSS_.BuildTree_FaceSplitting(Ax, Ay, px, py);
A.AFinv_HSS_.Construct_FaceSplitting(Ax, Ay);

mem_before = A.Storage();
rank_before = A.Rank();

A.AFinv_HSS_.Recompress(tol_hss);

mem_after = A.Storage();
rank_after = A.Rank();

compress_ratio = mem_after / mem_before;

fprintf("HSS compression: rank %d -> %d, stored bytes %.1e -> %.1e\n", ...
    rank_before, rank_after, mem_before, mem_after);
fprintf("compress ratio (after/before): %.1e\n", compress_ratio);

A.xy_inv_perm_ = zeros(size(A.xy_perm_));
A.xy_inv_perm_(A.xy_perm_) = (1 : size(A.xy_perm_, 1))';
A.omega_inv_perm_ = zeros(size(A.omega_perm_));
A.omega_inv_perm_(A.omega_perm_) = (1 : size(A.omega_perm_, 1))';

end
