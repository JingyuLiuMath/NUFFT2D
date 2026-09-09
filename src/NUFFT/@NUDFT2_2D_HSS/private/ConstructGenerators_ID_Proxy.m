function data = ConstructGenerators_ID_Proxy(A, xy, nx, ny, ranges, tol)
% Construct generators recursively using points in global leaf order.

local_nx = ranges(2) - ranges(1) + 1;
local_ny = ranges(4) - ranges(3) + 1;
data = struct;

if A.leaf_ == 1
    data.row_xy = xy(A.row_offset_ + (1 : A.row_size_), :);
    data.col_pos = TensorProduct2D((ranges(1) : ranges(2))', ...
        (ranges(3) : ranges(4))');
else
    children = cell(1, A.num_children_);
    x_child_size = local_nx;
    y_child_size = local_ny;
    if local_nx > 1
        x_child_size = [floor(local_nx / 2), ceil(local_nx / 2)];
    end
    if local_ny > 1
        y_child_size = [floor(local_ny / 2), ceil(local_ny / 2)];
    end
    num_children_x = size(x_child_size, 2);
    row_size = 0;
    col_size = 0;
    for i = 1 : A.num_children_
        ix = mod(i - 1, num_children_x) + 1;
        iy = floor((i - 1) / num_children_x) + 1;
        x_start = ranges(1) + (ix - 1) * x_child_size(1);
        y_start = ranges(3) + (iy - 1) * y_child_size(1);
        child_ranges = [x_start, x_start + x_child_size(ix) - 1, ...
            y_start, y_start + y_child_size(iy) - 1];
        children{i} = ConstructGenerators_ID_Proxy( ...
            A.children_{i}, xy, nx, ny, child_ranges, tol);
        row_size = row_size + A.children_{i}.row_rank_;
        col_size = col_size + A.children_{i}.col_rank_;
    end

    A.Bmat_ = cell(A.num_children_);
    for i = 1 : A.num_children_
        z_Ix = exp(-2 * pi * 1i * children{i}.row_xy(:, 1));
        z_Iy = exp(-2 * pi * 1i * children{i}.row_xy(:, 2));
        for j = 1 : A.num_children_
            if i ~= j
                w_Jx = exp(-2 * pi * 1i * children{j}.col_pos(:, 1) / nx);
                w_Jy = exp(-2 * pi * 1i * children{j}.col_pos(:, 2) / ny);
                A.Bmat_{i, j} = NUFFT2_Kernel(z_Ix, w_Jx, nx) ...
                    .* NUFFT2_Kernel(z_Iy, w_Jy, ny);
            end
        end
    end
    if A.level_ == 0
        return;
    end

    data.row_xy = zeros(row_size, 2);
    data.col_pos = zeros(col_size, 2);
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        data.row_xy(row_offset + (1 : A.children_{i}.row_rank_), :) = children{i}.row_xy;
        data.col_pos(col_offset + (1 : A.children_{i}.col_rank_), :) = children{i}.col_pos;
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end
end

z_Ix = exp(-2 * pi * 1i * data.row_xy(:, 1));
z_Iy = exp(-2 * pi * 1i * data.row_xy(:, 2));
w_Jx = exp(-2 * pi * 1i * data.col_pos(:, 1) / nx);
w_Jy = exp(-2 * pi * 1i * data.col_pos(:, 2) / ny);
if A.level_ == 0
    A.Amat_ = NUFFT2_Kernel(z_Ix, w_Jx, nx) .* NUFFT2_Kernel(z_Iy, w_Jy, ny);
    return;
end

% Use the same proxy layers and sample counts at every nonroot node.
local_n = max(local_nx, local_ny);
proxy_layer_size = 2 + log2(local_n);
rank_1d = ceil(1.5 * log(1 / tol) * log(2 * pi * (2 * local_n + 1)) * 2 / pi^2);
sampling_size_cross = ceil(2 * rank_1d * local_n);
sampling_size_diag = rank_1d^2;
proxy_surface = Proxy_Surface(nx, ny, ranges, proxy_layer_size, ...
    sampling_size_cross, sampling_size_diag);

% Construct row and column interpolation bases against the proxy surface.
A_I_proxy = NUFFT2_Kernel(z_Ix, proxy_surface(:, 1), nx) ...
    .* NUFFT2_Kernel(z_Iy, proxy_surface(:, 2), ny);
[row_sk, U, A.row_rank_] = LowRank_Row_ID(A_I_proxy, tol);
A_proxy_J = NUFFT2_Kernel(proxy_surface(:, 1), w_Jx, nx) ...
    .* NUFFT2_Kernel(proxy_surface(:, 2), w_Jy, ny);
[col_sk, V, A.col_rank_] = LowRank_ID(A_proxy_J, tol);

if A.row_offset_ == 0 && A.col_offset_ == 0
    fprintf("    \n");
    fprintf("    level: %d\n", A.level_);
    fprintf("    local_n: %d\n", local_n);
    fprintf("    row size: %d, col size: %d\n", size(data.row_xy, 1), size(data.col_pos, 1));
    fprintf("    rank_1d: %d\n", rank_1d);
    fprintf("    sampling_size_cross: %d\n", sampling_size_cross);
    fprintf("    sampling_size_diag: %d\n", sampling_size_diag);
    fprintf("    row rank: %d, col rank: %d\n", A.row_rank_, A.col_rank_);
end

if A.leaf_ == 1
    A.Umat_ = U;
    A.Vmat_ = V;
    A.Amat_ = NUFFT2_Kernel(z_Ix, w_Jx, nx) .* NUFFT2_Kernel(z_Iy, w_Jy, ny);
else
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        A.children_{i}.Rmat_ = U(row_offset + (1 : A.children_{i}.row_rank_), :);
        A.children_{i}.Wmat_ = V(col_offset + (1 : A.children_{i}.col_rank_), :);
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end
end
data.row_xy = data.row_xy(row_sk, :);
data.col_pos = data.col_pos(col_sk, :);

end
