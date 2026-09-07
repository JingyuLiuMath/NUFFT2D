function [Q, k] = ColBasis(B, k, tol)
% ColBasis

arguments (Input)
    B (:, :) double;
    k (1, 1) double;
    tol (1, 1) double;
end

arguments (Output)
    Q (:, :) double;
    k (1, 1) double;
end

[Q, ~, ~, k1] = MyQRSketch(B, tol);
k = min(k, k1);
Q = Q(:, 1 : k);

end