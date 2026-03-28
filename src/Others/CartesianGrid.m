function x = CartesianGrid(n, alpha)

arguments (Input)
    n (1, 1) double;
    alpha (1, 1) double = 2;
end

arguments (Output)
    x (:, 2) double;
end

m = alpha * n;

x = (0 : (m - 1))' / m;
x = TensorProduct2D(x, x);

end