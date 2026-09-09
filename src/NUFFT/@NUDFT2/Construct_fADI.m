function Construct_fADI(A, x, n_leaf, tol)

arguments (Input)
    A NUDFT2;
    x (:, 1) double;
    n_leaf (1, 1) double;
    tol (1, 1) double;
end

A.Factor_ = struct();

A.M_ = size(x, 1);
M = A.M_;
N = A.N_;

A.AFinv_HSS_ = NUDFT2_HSS(N);
eta = ceil(M / N);
A.x_perm_ = A.AFinv_HSS_.BuildTree(x, n_leaf, eta);
A.x_inv_perm_ = zeros(size(A.x_perm_));
A.x_inv_perm_(A.x_perm_) = (1 : size(A.x_perm_, 1))';

A.AFinv_HSS_.Construct_fADI(tol);

end
