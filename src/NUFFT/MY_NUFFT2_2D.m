function f = MY_NUFFT2_2D(u, xy, nx, ny, method, acc)
% MY_NUFFT2_2D

arguments (Input)
    u (:, :) double;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    method string = "finufft";
    acc (1, 1) double = 1e-15;
end


arguments (Output)
    f (:, :) double;
end

switch method
    case "direct"
        f = NUDFT2_2D_Matrix(xy, nx, ny) * u;
    case "finufft"
        u = reshape(u, nx, ny ,[]);
        f = finufft2d2(...
            2 * pi * xy(:, 1), 2 * pi * xy(:, 2), ...
            -1, acc, u);
        d_offset = ...
            exp(-2 * pi * 1i * floor(nx / 2) * xy(:, 1)) ...
            .* exp(-2 * pi * 1i * floor(ny / 2) * xy(:, 2));
        f = d_offset .* f;
end

end