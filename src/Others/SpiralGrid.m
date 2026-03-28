function x = SpiralGrid(M, circle_num)

arguments (Input)
    M (1, 1) double;
    circle_num = 20;
end

arguments (Output)
    x (:, 2) double;
end

c = 0.5 + 1i * 0.5;
t = linspace(0, circle_num, M)';
z = c + sqrt(2) / 2 * t / circle_num .* exp(2 * pi * 1i * t);
x = [real(z), imag(z)];
x = x(x(:, 1) >= 0, :);
x = x(x(:, 1) < 1, :);
x = x(x(:, 2) >= 0, :);
x = x(x(:, 2) < 1, :);

end