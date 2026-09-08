function proxy_surface = Proxy_Surface(nx, ny, ranges, proxy_layer_size, ...
    sampling_size_cross, sampling_size_diag)
% Sample the four strips and four corners surrounding a column box.

hx = 1 / nx;
hy = 1 / ny;
x_start = ranges(1) / nx - hx / 2;
x_end = ranges(2) / nx + hx / 2;
y_start = ranges(3) / ny - hy / 2;
y_end = ranges(4) / ny + hy / 2;
dx = proxy_layer_size * hx;
dy = proxy_layer_size * hy;
bounds = [x_start, x_end, y_start - dy, y_start; ...
    x_start, x_end, y_end, y_end + dy; ...
    x_start - dx, x_start, y_start, y_end; ...
    x_end, x_end + dx, y_start, y_end; ...
    x_start - dx, x_start, y_start - dy, y_start; ...
    x_start - dx, x_start, y_end, y_end + dy; ...
    x_end, x_end + dx, y_start - dy, y_start; ...
    x_end, x_end + dx, y_end, y_end + dy];
proxy_surface = zeros(4 * sampling_size_cross + 4 * sampling_size_diag, 2);
row_offset = 0;
for i = 1 : 8
    if i <= 4
        sampling_size = sampling_size_cross;
    else
        sampling_size = sampling_size_diag;
    end
    proxy_surface(row_offset + (1 : sampling_size), :) = RandRectangular( ...
        sampling_size, bounds(i, 1), bounds(i, 2), bounds(i, 3), bounds(i, 4));
    row_offset = row_offset + sampling_size;
end
proxy_surface = exp(-2 * pi * 1i * proxy_surface);

end
