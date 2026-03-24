function p = BuildTree(A, row_x, min_points, eta, store_points)
% BuildTree

% Jingyu Liu, November 24, 2024.

arguments (Input)
    A NUDFT2_HSS;
    row_x (:, 1) double;
    min_points (1, 1) double;
    eta (1, 1) double;
    store_points (1, 1) double = 0;
end

arguments (Output)
    p (:, 1) double;  % Permutation such that x(p) follows the leaf order.
end

x_length = size(row_x, 1);
A.row_size_ = x_length;

if A.level_ == 0
    A.row_global_size_ = x_length;
end

if A.col_size_ <= min_points || A.row_size_ <= eta * min_points
    A.leaf_ = 1;
    A.max_level_ = A.level_;
    p = (1 : x_length)';
    if store_points == 1
        A.row_x_ = row_x;
        A.col_pos_ = (A.pos_start_ : A.pos_end_)';
    end
else
    % Partition.
    A.num_children_ = 2;
    A.children_ = cell(1, A.num_children_);
    c1_size = floor(A.col_size_ / 2);
    child_size = [c1_size, A.col_size_ - c1_size];

    p = zeros(x_length, 1);
    I = cell(1, A.num_children_);
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        current_col_size = child_size(i);
        A.children_{i} = NUDFT2_HSS(...
            A.col_global_size_, ...
            A.pos_start_ + col_offset, ...
            A.pos_start_ + col_offset + current_col_size - 1, ...
            A.level_ + 1, A.level_order_ * 2 + i - 1, ...
            A.node_order_ * 2 + i);

        I{i} = FindID_ExtendArc(...
            A.col_global_size_, ...
            A.children_{i}.pos_start_, ...
            A.children_{i}.pos_end_, ...
            row_x);
        current_row_size = size(I{i}, 1);
        p((row_offset + 1) : (row_offset + current_row_size)) = I{i};

        row_offset = row_offset + current_row_size;
        col_offset = col_offset + current_col_size;
    end

    % Recursion.
    row_offset = 0;
    for i = 1 : A.num_children_
        A.children_{i}.row_global_size_ = A.row_global_size_;
        pi = A.children_{i}.BuildTree(...
            row_x(I{i}), min_points, eta, store_points);
        current_row_size = length(I{i});
        c_pi = p((row_offset + 1) : (row_offset + current_row_size));
        p((row_offset + 1) : (row_offset + current_row_size)) = c_pi(pi);
        row_offset = row_offset + current_row_size;
        A.max_level_ = max(A.max_level_, A.children_{i}.max_level_);
    end
end

end