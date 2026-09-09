function [p, q] = BuildTree_FaceSplitting(A, Ax, Ay, px, py, data)

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    px (:, 1) double;
    py (:, 1) double;
    data struct = struct();
end

arguments (Output)
    p (:, 1) double;
    q (:, 1) double;
end

p = zeros(0, 1);
q = zeros(0, 1);

if A.level_ == 0
    % Group rows by their pair of 1D leaves.
    M = size(px, 1);
    x_inv_perm = zeros(M, 1);
    y_inv_perm = zeros(M, 1);
    x_inv_perm(px) = (1 : M)';
    y_inv_perm(py) = (1 : M)';
    [row_leaf_x, col_leaf_x, num_leaves_x] = LeafRows(Ax);
    [row_leaf_y, col_leaf_y, num_leaves_y] = LeafRows(Ay);
    leaf_pair = row_leaf_x(x_inv_perm) ...
        + num_leaves_x * (row_leaf_y(y_inv_perm) - 1);
    num_leaf_pairs = num_leaves_x * num_leaves_y;
    row_counts = accumarray(leaf_pair, ones(M, 1), ...
        [num_leaf_pairs, 1]);
    row_start = cumsum([0; row_counts]);
    row_offset = row_start(1 : end - 1);

    data.p = zeros(M, 1);
    data.row_xy = zeros(M, 2);
    for i = 1 : M
        leaf = leaf_pair(i);
        row_offset(leaf) = row_offset(leaf) + 1;
        data.p(row_offset(leaf)) = i;
        data.row_xy(row_offset(leaf), :) ...
            = [x_inv_perm(i), y_inv_perm(i)];
    end
    data.row_start = row_start;
    data.col_leaf_x = col_leaf_x;
    data.col_leaf_y = col_leaf_y;
    data.num_leaves_x = num_leaves_x;
    row_counts = reshape(row_counts, num_leaves_x, num_leaves_y);
    data.row_prefix = zeros(num_leaves_x + 1, num_leaves_y + 1);
    data.row_prefix(2 : end, 2 : end) ...
        = cumsum(cumsum(row_counts, 1), 2);
    data.root = A;

    A.row_global_size_ = M;
    A.row_xy_ = zeros(M, 2);
    A.x_col_pos_ = zeros(A.col_size_, 1);
end

% Current node.
A.max_level_ = A.level_;
x_leaf_begin = data.col_leaf_x(Ax.col_offset_ + 1);
x_leaf_end = data.col_leaf_x(Ax.col_offset_ + Ax.col_size_);
y_leaf_begin = data.col_leaf_y(Ay.col_offset_ + 1);
y_leaf_end = data.col_leaf_y(Ay.col_offset_ + Ay.col_size_);
A.row_size_ = data.row_prefix(x_leaf_end + 1, y_leaf_end + 1) ...
    - data.row_prefix(x_leaf_begin, y_leaf_end + 1) ...
    - data.row_prefix(x_leaf_end + 1, y_leaf_begin) ...
    + data.row_prefix(x_leaf_begin, y_leaf_begin);

if Ax.leaf_ == 1 && Ay.leaf_ == 1
    % Leaf node.
    A.leaf_ = 1;
    leaf = x_leaf_begin + data.num_leaves_x * (y_leaf_begin - 1);
    I = data.row_start(leaf) + (1 : A.row_size_);
    row_xy = data.row_xy(I, :);
    row_xy(:, 1) = row_xy(:, 1) - Ax.row_offset_;
    row_xy(:, 2) = row_xy(:, 2) - Ay.row_offset_;
    q_leaf = reshape((Ax.pos_start_ : Ax.pos_end_)' ...
        + A.nx_ * (Ay.pos_start_ : Ay.pos_end_) + 1, [], 1);

    if A.level_ == 0
        p = data.p(I);
        q = q_leaf;
        A.row_xy_ = row_xy;
        A.x_col_pos_ = [];
    else
        row_index = A.row_offset_ + (1 : A.row_size_);
        data.root.row_xy_(row_index, 1) = data.p(I);
        A.row_xy_ = row_xy;
        col_index = A.col_offset_ + (1 : A.col_size_);
        data.root.x_col_pos_(col_index) = q_leaf;
    end
    return;
end


% Partition.
num_children_x = max(Ax.num_children_, 1);
num_children_y = max(Ay.num_children_, 1);
A.num_children_ = num_children_x * num_children_y;
A.children_ = cell(1, A.num_children_);
row_offset = 0;
col_offset = 0;
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;

    if Ax.leaf_ == 1
        x_pos_start = Ax.pos_start_;
        x_pos_end = Ax.pos_end_;
    else
        x_pos_start = Ax.children_{ix}.pos_start_;
        x_pos_end = Ax.children_{ix}.pos_end_;
    end
    if Ay.leaf_ == 1
        y_pos_start = Ay.pos_start_;
        y_pos_end = Ay.pos_end_;
    else
        y_pos_start = Ay.children_{iy}.pos_start_;
        y_pos_end = Ay.children_{iy}.pos_end_;
    end

    A.children_{i} = NUDFT2_2D_HSS(...
        A.nx_, A.ny_, ...
        x_pos_start, x_pos_end, y_pos_start, y_pos_end, ...
        A.level_ + 1, A.row_offset_ + row_offset, ...
        A.col_offset_ + col_offset);
    A.children_{i}.row_global_size_ = A.row_global_size_;

    if Ax.leaf_ == 1
        A.children_{i}.BuildTree_FaceSplitting(...
            Ax, Ay.children_{iy}, px, py, data);
    elseif Ay.leaf_ == 1
        A.children_{i}.BuildTree_FaceSplitting(...
            Ax.children_{ix}, Ay, px, py, data);
    else
        A.children_{i}.BuildTree_FaceSplitting(...
            Ax.children_{ix}, Ay.children_{iy}, px, py, data);
    end

    row_offset = row_offset + A.children_{i}.row_size_;
    col_offset = col_offset + A.children_{i}.col_size_;
    A.max_level_ = max(A.max_level_, A.children_{i}.max_level_);
end

if A.level_ == 0
    p = A.row_xy_(:, 1);
    q = A.x_col_pos_;
    A.row_xy_ = [];
    A.x_col_pos_ = [];
end

end
