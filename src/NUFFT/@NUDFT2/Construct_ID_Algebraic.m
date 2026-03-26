function Construct_ID_Algebraic(A, x, min_points, rank_or_tol)

arguments (Input)
    A NUDFT2;
    x (:, 1) double;
    min_points (1, 1) double;
    rank_or_tol (1, 1) double;
end

A.M_ = size(x, 1);
M = A.M_;
N = A.N_;

A.AFinv_HSS_ = NUDFT2_HSS(N);
eta = ceil(M / N);
A.x_perm_ = A.AFinv_HSS_.BuildTree(x, min_points, eta, 0);
[~, A.x_inv_perm_] = sort(A.x_perm_, "ascend");

H = NUDFT2_Matrix(x, N);
H = ifft(H, [], 2);
H = H(A.x_perm_, :);

A.AFinv_HSS_.Construct_ID_Algebraic(H, rank_or_tol);

end