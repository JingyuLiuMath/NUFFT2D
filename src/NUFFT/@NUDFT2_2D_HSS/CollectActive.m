function [row_active, col_active] = CollectActive(A, level)

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
end

arguments (Output)
    row_active (:, 1) double;
    col_active (:, 1) double;
end

if A.level_ == level || A.leaf_ == 1
    row_active = A.row_ind_;
    col_active = A.col_ind_;
else
    row_sk_offset = A.row_sk_offset_;
    col_sk_offset = A.col_sk_offset_;
    row_active = [];
    col_active = [];
    for i = 1 : A.num_children_
        A.children_{i}.row_sk_offset_ = row_sk_offset;
        A.children_{i}.col_sk_offset_ = col_sk_offset;
        [tmp_row, tmp_col] = A.children_{i}.CollectActive(level);
        row_active = [row_active; tmp_row];
        col_active = [col_active; tmp_col];
        row_sk_offset = row_sk_offset + size(tmp_row, 1);
        col_sk_offset = col_sk_offset + size(tmp_col, 1);
    end
end

A.row_sk_size_= size(row_active, 1);
A.col_sk_size_= size(col_active, 1);

end