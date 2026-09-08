function Construct_fADI(A, x, n_leaf, tol)

arguments (Input)
    A NUDFT2;
    x (:, 1) double;
    n_leaf (1, 1) double;
    tol (1, 1) double;
end

A.Factor_ = struct();
A.M_ = size(x, 1);
N = A.N_;

A.AFinv_HSS_ = NUDFT2_HSS(N);
A.x_perm_ = A.AFinv_HSS_.BuildTree(x, n_leaf);
A.AFinv_HSS_.Construct_fADI(x(A.x_perm_), tol);
A.x_inv_perm_ = zeros(size(A.x_perm_));
A.x_inv_perm_(A.x_perm_) = (1 : size(A.x_perm_, 1))';

end
