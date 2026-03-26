function [row_inact, col_inact] = Deactivate(A, level, row_inact, col_inact)

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
    row_inact (:, 1) double;
    col_inact (:, 1) double;
end

arguments (Output)
    row_inact (:, 1) double;
    col_inact (:, 1) double;
end

if A.level_ == level
    row_inact = [row_inact; A.row_re_];
    col_inact = [col_inact; A.col_re_];
    A.row_re_ = [];
    A.col_re_ = [];
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        [row_inact, col_inact] = A.children_{i}.Deactivate(level, row_inact, col_inact);
    end
end

end