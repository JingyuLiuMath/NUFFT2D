function x = SpiralGrid(N, num_curve, circle_num)

arguments (Input)
    N (1, 1) double;
    num_curve = ceil(log2(N) / 2);
    circle_num = 20;
end

arguments (Output)
    x (:, 2) double;
end

c = 0.5 + 1i * 0.5;
base_radius = sqrt(2) / 2;

x = [];

for k = 1:num_curve
    angle_offset = 2 * pi * (k - 1) / num_curve;

    t = linspace(0, circle_num, N)';

    z = c + base_radius * t / circle_num .* exp(2 * pi * 1i * t + 1i * angle_offset);

    curve_points = [real(z), imag(z)];

    curve_points = curve_points(curve_points(:, 1) >= 0, :);
    curve_points = curve_points(curve_points(:, 1) < 1, :);
    curve_points = curve_points(curve_points(:, 2) >= 0, :);
    curve_points = curve_points(curve_points(:, 2) < 1, :);

    x = [x; curve_points];
end

end