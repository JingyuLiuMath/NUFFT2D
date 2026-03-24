function [row_active, col_active] = CollectActive(A, level)

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
end

arguments (Output)
    row_active (:, 1) double;
    col_active (:, 1) double;
end

if A.level_ == level
    row_active = A.row_ind_;
    col_active = A.col_ind_;
else
    row_active = [];
    col_active = [];
    for i = 1 : A.num_children_
        [tmp_row, tmp_col] = A.children_{i}.CollectActive(level);
        row_active = [row_active; tmp_row];
        col_active = [col_active; tmp_col];
    end
end

end