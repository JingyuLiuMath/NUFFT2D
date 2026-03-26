function f = MY_NUFFT1_2D(u, omega, mx, my, method, acc)
% MY_NUFFT1_2D

% f[k] = \sum_{jx = 0}^{mx - 1} \sum_{jy = 0}^{my - 1} ...
%     exp(-2 * pi * 1i * (jx / mx) * omega[k, 1]) ...
%     * exp(-2 * pi * 1i * (jy / my) * omega[k, 2]) ...
%     * u[jx, jy]

arguments (Input)
    u (:, :) double;
    omega (:, 2) double;
    mx (1, 1) double;
    my (1, 1) double;
    method string = "finufft";
    acc (1, 1) double = 1e-15;
end


arguments (Output)
    f (:, :) double;
end

switch method
    case "direct"
        f = NUDFT1_2D_Matrix(omega, mx, my) * u;
    case "finufft"
        d_offset = ...
            exp(-2 * pi * 1i * floor(mx / 2) * omega(:, 1) / mx) ...
            .* exp(-2 * pi * 1i * floor(my / 2) * omega(:, 2) / my);
        u = d_offset .* u;

        f = finufft2d1(...
            2 * pi * omega(:, 1) / mx,...
            2 * pi * omega(:, 2) / my, ...
            u,-1,acc,mx,my);

        f = reshape(f, mx * my, []);
end

end