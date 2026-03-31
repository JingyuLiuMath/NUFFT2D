function xy = RandRectangular(M, a, b, c, d)

arguments (Input)
    M (1, 1) double;
    a (1, 1) double;
    b (1, 1) double;
    c (1, 1) double;
    d (1, 1) double;
end

arguments(Output)
    xy (:, 2) double;
end

xy = rand(M, 2);

xy(:, 1) = xy(:, 1) * (b - a) + a;
xy(:, 2) = xy(:, 2) * (d - c) + c;

end