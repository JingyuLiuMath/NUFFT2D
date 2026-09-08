function Construct_FaceSplitting(A, xy, n_leaf, tol_hss)

arguments (Input)
    A NUDFT2_2D;
    xy (:, 2) double;
    n_leaf (1, 1) double;
    tol_hss (1, 1) double;
end

A.Factor_ = struct();
Ax = NUDFT2(A.nx_);
Ax.Construct_fADI(xy(:, 1), n_leaf, tol_hss / 10);
Ay = NUDFT2(A.ny_);
Ay.Construct_fADI(xy(:, 2), n_leaf, tol_hss / 10);

A.AFinv_HSS_ = NUDFT2_2D_HSS(A.nx_, A.ny_);
[p, q] = A.AFinv_HSS_.Construct_FaceSplitting( ...
    Ax.AFinv_HSS_, Ay.AFinv_HSS_, Ax.x_perm_, Ay.x_perm_);
mem_before = A.Storage();
rank_before = A.Rank();
A.AFinv_HSS_.Recompress(tol_hss);
mem_after = A.Storage();
compress_ratio = mem_after / mem_before;
fprintf("HSS compression: rank %d -> %d, stored bytes %.1e -> %.1e\n", ...
    rank_before, A.Rank(), mem_before, mem_after);
fprintf("compress ratio (after/before): %.1e\n", compress_ratio);
A.M_ = size(xy, 1);
A.xy_perm_ = p;
A.omega_perm_ = q;
A.xy_inv_perm_ = zeros(size(p));
A.xy_inv_perm_(p) = (1 : numel(p))';
A.omega_inv_perm_ = zeros(size(q));
A.omega_inv_perm_(q) = (1 : numel(q))';

end
