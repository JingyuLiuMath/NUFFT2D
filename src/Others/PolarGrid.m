function x = PolarGrid(n, alpha_radius, alpha_t)

arguments (Input)
    n (1, 1) double;
    alpha_radius (1, 1) double = 1;
    alpha_t (1, 1) double = 4;
end

arguments (Output)
    x (:, 2) double;
end

n_radius = ceil(alpha_radius * n);
n_t = ceil(alpha_t * n);

c = 0.5 + 1i * 0.5;
r = sqrt(2) / 2 * (1 : (n_radius - 1))' / n_radius;
t = (0 : (n_t - 1))' / (n_t);
z = c + r .* exp(2 * pi * 1i * t');
z = z(:);
x = [real(z), imag(z)];
x = x(x(:, 1) >= 0, :);
x = x(x(:, 1) < 1, :);
x = x(x(:, 2) >= 0, :);
x = x(x(:, 2) < 1, :);
x = [x; real(c), imag(c)];

end