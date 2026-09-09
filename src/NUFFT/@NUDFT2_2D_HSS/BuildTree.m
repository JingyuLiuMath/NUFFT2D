function [p, q] = BuildTree(A, xy, nx, ny, N_leaf, ranges)
% BuildTree Build the tree and return row and column permutations.

arguments (Input)
    A NUDFT2_2D_HSS;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    N_leaf (1, 1) double;
    ranges (1, 4) double = [0, nx - 1, 0, ny - 1];
end

arguments (Output)
    p (:, 1) double;
    q (:, 1) double;
end

A.row_size_ = size(xy, 1);
if A.level_ == 0
    A.row_global_size_ = A.row_size_;
end
eta = ceil(A.row_global_size_ / A.col_global_size_);
A.max_level_ = A.level_;
local_nx = ranges(2) - ranges(1) + 1;
local_ny = ranges(4) - ranges(3) + 1;

if A.col_size_ <= N_leaf || A.row_size_ <= eta * N_leaf
    A.leaf_ = 1;
    p = (1 : A.row_size_)';
    col_pos = TensorProduct2D((ranges(1) : ranges(2))', (ranges(3) : ranges(4))');
    q = sub2ind([nx, ny], col_pos(:, 1) + 1, col_pos(:, 2) + 1);
else
    x_child_size = local_nx;
    y_child_size = local_ny;
    if local_nx > 1
        x_child_size = [floor(local_nx / 2), ceil(local_nx / 2)];
    end
    if local_ny > 1
        y_child_size = [floor(local_ny / 2), ceil(local_ny / 2)];
    end
    num_children_x = size(x_child_size, 2);
    A.num_children_ = num_children_x * size(y_child_size, 2);
    A.children_ = cell(1, A.num_children_);
    p = zeros(A.row_size_, 1);
    q = zeros(A.col_size_, 1);
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        ix = mod(i - 1, num_children_x) + 1;
        iy = floor((i - 1) / num_children_x) + 1;
        x_start = ranges(1) + (ix - 1) * x_child_size(1);
        y_start = ranges(3) + (iy - 1) * y_child_size(1);
        child_ranges = [x_start, x_start + x_child_size(ix) - 1, ...
            y_start, y_start + y_child_size(iy) - 1];
        A.children_{i} = NUDFT2_2D_HSS(nx, ny);
        A.children_{i}.col_size_ = x_child_size(ix) * y_child_size(iy);
        A.children_{i}.level_ = A.level_ + 1;
        A.children_{i}.row_global_size_ = A.row_global_size_;
        A.children_{i}.row_offset_ = A.row_offset_ + row_offset;
        A.children_{i}.col_offset_ = A.col_offset_ + col_offset;
        Ix = FindID_ExtendArc(nx, child_ranges(1), child_ranges(2), xy(:, 1));
        Iy = FindID_ExtendArc(ny, child_ranges(3), child_ranges(4), xy(:, 2));
        I = intersect(Ix, Iy, "sorted");
        [p_i, q_i] = A.children_{i}.BuildTree(xy(I, :), nx, ny, N_leaf, child_ranges);
        p(row_offset + (1 : size(I, 1))) = I(p_i);
        q(col_offset + (1 : A.children_{i}.col_size_)) = q_i;
        row_offset = row_offset + size(I, 1);
        col_offset = col_offset + A.children_{i}.col_size_;
        A.max_level_ = max(A.max_level_, A.children_{i}.max_level_);
    end
end

end
