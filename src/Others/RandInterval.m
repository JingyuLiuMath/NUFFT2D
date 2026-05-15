function x = RandInterval(M, a, b)

arguments (Input)
    M (1, 1) double;
    a (1, 1) double;
    b (1, 1) double;
end

arguments(Output)
    x (:, 1) double;
end

x = rand(M, 1);

x = x * (b - a) + a;

end