function Construct_BlackBox(A, x, min_points, target_rank, tol)
% BlackBoxConstruct

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A NUDFT2;
    x (:, 1) double;
    min_points (1, 1) double;
    target_rank (1, 1) double;
    tol (1, 1) double;
end

A.M_ = size(x, 1);
N = A.N_;
M = A.M_;

A.AFinv_HSS_ = NUDFT2_HSS(N);
eta = ceil(M / N);
A.x_perm_ = A.AFinv_HSS_.BuildTree(x, min_points, eta, 0);
[~, A.x_inv_perm_] = sort(A.x_perm_, "ascend");
[row_leaf_size, col_leaf_size] = A.AFinv_HSS_.MaxLeafSize();

x_p = x(A.x_perm_);
op_A = @(v) MY_NUFFT2(ifft(v), x_p, N);
op_Astar = @(v) fft(MY_NUFFT1(v, -N * x_p, N)) / N;
A.AFinv_HSS_.BlackBoxConstruct(op_A, op_Astar, ...
    row_leaf_size, col_leaf_size, target_rank, tol);

end