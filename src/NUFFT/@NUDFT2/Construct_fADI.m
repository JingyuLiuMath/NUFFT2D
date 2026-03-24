function Construct_fADI(A, x, min_points, tol)
% Construct

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A NUDFT2;
    x (:, 1) double;
    min_points (1, 1) double;
    tol (1, 1) double;
end

A.M_ = size(x, 1);
M = A.M_;
N = A.N_;

A.AFinv_HSS_ = NUDFT2_HSS(N);
eta = ceil(M / N);
A.x_perm_ = A.AFinv_HSS_.BuildTree(x, min_points, eta, 1);
[~, A.x_inv_perm_] = sort(A.x_perm_, "ascend");
A.AFinv_HSS_.Construct_fADI(tol);

end