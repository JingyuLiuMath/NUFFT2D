function [p, q] = BuildTree(A, row_xy, min_points, eta, store_points)
% BuildTree

% Jingyu Liu, January 27, 2025.

arguments (Input)
    A NUDFT2_2D_HSS;
    row_xy (:, 2) double;
    min_points (1, 1) double;
    eta (1, 1) double;
    store_points (1, 1) double = 0;
end

arguments (Output)
    p (:, 1) double;
    q (:, 1) double;
end

xy_length = size(row_xy, 1);
A.row_size_ = xy_length;

if A.level_ == 0
    A.row_global_size_ = xy_length;
end

if A.col_size_ <= min_points || A.row_size_ <= eta * min_points
    A.leaf_ = 1;
    A.max_level_ = A.level_;

    x_col_pos = (A.x_pos_start_ : A.x_pos_end_)';
    y_col_pos = (A.y_pos_start_ : A.y_pos_end_)';
    xy_col_pos  = TensorProduct2D(x_col_pos, y_col_pos);

    p = (1 : xy_length)';
    q = sub2ind([A.nx_, A.ny_], ...
        xy_col_pos(:, 1) + 1, xy_col_pos(:, 2) + 1);

    if store_points == 1
        A.row_xy_ = row_xy;
        A.col_pos_ = xy_col_pos;
    end
else
    % Partition.
    A.num_children_ = 4;
    A.children_ = cell(1, A.num_children_);
    x_c1_size = floor(A.x_col_size_ / 2);
    x_child_size = [x_c1_size, A.x_col_size_ - x_c1_size];
    y_c1_size = floor(A.y_col_size_ / 2);
    y_child_size = [y_c1_size, A.y_col_size_ - y_c1_size];

    p = zeros(xy_length, 1);
    q = zeros(A.col_size_, 1);
    it = 0;
    row_offset = 0;
    y_col_offset = 0;
    I = cell(1, A.num_children_);
    for ity = 1 : 2
        curr_y_col_size = y_child_size(ity);
        x_col_offset = 0;
        for itx = 1 : 2
            it = it + 1;
            curr_x_col_size = x_child_size(itx);
            A.children_{it} = NUDFT2_2D_HSS(A.nx_, A.ny_, ...
                A.x_pos_start_ + x_col_offset, ...
                A.x_pos_start_ + x_col_offset + curr_x_col_size - 1, ...
                A.y_pos_start_ + y_col_offset, ...
                A.y_pos_start_ + y_col_offset + curr_y_col_size - 1, ...
                A.level_ + 1, ...
                A.level_order_ * 4 + it - 1, ...
                A.node_order_ * 4 + it);

            Ix = FindID_ExtendArc(...
                A.nx_, ...
                A.children_{it}.x_pos_start_, ...
                A.children_{it}.x_pos_end_, ...
                row_xy(:, 1));
            Iy = FindID_ExtendArc(...
                A.ny_, ...
                A.children_{it}.y_pos_start_, ...
                A.children_{it}.y_pos_end_, ...
                row_xy(:, 2));
            I{it} = intersect(Ix, Iy, "sorted");
            curr_row_size = size(I{it}, 1);
            p((row_offset + 1) : (row_offset + curr_row_size)) = I{it};

            row_offset = row_offset + curr_row_size;
            x_col_offset = x_col_offset + curr_x_col_size;
        end
        y_col_offset = y_col_offset + curr_y_col_size;
    end

    % Recursion.
    row_offset = 0;
    col_offset = 0;
    for it = 1 : A.num_children_
        A.children_{it}.row_global_size_ = A.row_global_size_;
        [pi, qi] = A.children_{it}.BuildTree(...
            row_xy(I{it}, :), min_points, eta, store_points);
        curr_row_size = length(I{it});
        curr_col_size = A.children_{it}.col_size_;
        c_pi = p((row_offset + 1) : (row_offset + curr_row_size));
        p((row_offset + 1) : (row_offset + curr_row_size)) = c_pi(pi);
        q((col_offset + 1) : (col_offset + curr_col_size)) = qi;
        A.max_level_ = max(A.max_level_, A.children_{it}.max_level_);

        row_offset = row_offset + curr_row_size;
        col_offset = col_offset + curr_col_size;
    end
end

end