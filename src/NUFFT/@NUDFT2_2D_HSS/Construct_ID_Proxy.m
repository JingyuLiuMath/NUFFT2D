function [p, q] = Construct_ID_Proxy(A, xy, nx, ny, N_leaf, tol)
% Construct_ID_Proxy Build a fresh 2D HSS root using proxy interpolation.
% Point coordinates and skeletons are local to the construction.

arguments (Input)
    A NUDFT2_2D_HSS;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    N_leaf (1, 1) double;
    tol (1, 1) double;
end

arguments (Output)
    p (:, 1) double;
    q (:, 1) double;
end

A.row_global_size_ = size(xy, 1);
eta = ceil(A.row_global_size_ / A.col_global_size_);
[p, q] = ConstructGenerators_ID_Proxy(A, xy, nx, ny, ...
    [0, nx - 1, 0, ny - 1], N_leaf, eta, tol);

end
