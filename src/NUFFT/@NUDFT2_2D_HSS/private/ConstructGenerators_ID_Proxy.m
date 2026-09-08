function [p, q, data] = ConstructGenerators_ID_Proxy(A, row_xy, nx, ny, ...
    ranges, N_leaf, eta, tol)
% Build the quadtree and pass interpolation skeletons to the parent.

A.row_size_ = size(row_xy, 1);
A.max_level_ = A.level_;
local_nx = ranges(2) - ranges(1) + 1;
local_ny = ranges(4) - ranges(3) + 1;
data = struct;

if A.col_size_ <= N_leaf || A.row_size_ <= eta * N_leaf
    A.leaf_ = 1;
    p = (1 : A.row_size_)';
    data.row_xy = row_xy;
    data.col_pos = TensorProduct2D((ranges(1) : ranges(2))', ...
        (ranges(3) : ranges(4))');
    q = sub2ind([nx, ny], data.col_pos(:, 1) + 1, data.col_pos(:, 2) + 1);
else
    A.num_children_ = 4;
    A.children_ = cell(1, A.num_children_);
    children = cell(1, A.num_children_);
    x_first_size = floor(local_nx / 2);
    y_first_size = floor(local_ny / 2);
    x_child_size = [x_first_size, local_nx - x_first_size];
    y_child_size = [y_first_size, local_ny - y_first_size];
    p = zeros(A.row_size_, 1);
    q = zeros(A.col_size_, 1);
    row_offset = 0;
    col_offset = 0;
    row_size = 0;
    col_size = 0;
    for i = 1 : A.num_children_
        ix = mod(i - 1, 2) + 1;
        iy = floor((i - 1) / 2) + 1;
        x_start = ranges(1) + (ix - 1) * x_first_size;
        y_start = ranges(3) + (iy - 1) * y_first_size;
        child_ranges = [x_start, x_start + x_child_size(ix) - 1, ...
            y_start, y_start + y_child_size(iy) - 1];
        A.children_{i} = NUDFT2_2D_HSS(nx, ny);
        A.children_{i}.col_size_ = x_child_size(ix) * y_child_size(iy);
        A.children_{i}.level_ = A.level_ + 1;
        A.children_{i}.row_global_size_ = A.row_global_size_;
        A.children_{i}.row_offset_ = A.row_offset_ + row_offset;
        A.children_{i}.col_offset_ = A.col_offset_ + col_offset;
        Ix = FindID_ExtendArc(nx, child_ranges(1), child_ranges(2), row_xy(:, 1));
        Iy = FindID_ExtendArc(ny, child_ranges(3), child_ranges(4), row_xy(:, 2));
        I = intersect(Ix, Iy, "sorted");
        [p_i, q_i, children{i}] = ConstructGenerators_ID_Proxy( ...
            A.children_{i}, row_xy(I, :), nx, ny, child_ranges, N_leaf, eta, tol);
        p(row_offset + (1 : numel(I))) = I(p_i);
        q(col_offset + (1 : A.children_{i}.col_size_)) = q_i;
        row_offset = row_offset + numel(I);
        col_offset = col_offset + A.children_{i}.col_size_;
        row_size = row_size + A.children_{i}.row_rank_;
        col_size = col_size + A.children_{i}.col_rank_;
        A.max_level_ = max(A.max_level_, A.children_{i}.max_level_);
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
