function AssignHSSRowColumn(A, row_ind, col_ind)
% AssignHSSRowColumn

% Jingyu Liu, November 18, 2024.

arguments (Input)
    A NUDFT2_HSS;
    row_ind (:, 1) double;
    col_ind (:, 1) double;
end

A.self_row_ind_ = row_ind;
A.self_col_ind_ = col_ind;

if A.leaf_ == 1
    A.row_ind_ = row_ind;
    A.col_ind_ = col_ind;
end

if A.leaf_ == 0
    % Recursion.
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        curr_row_size = A.children_{i}.row_size_;
        curr_col_size = A.children_{i}.col_size_;
        A.children_{i}.AssignHSSRowColumn(...
            row_ind((row_offset + 1) : (row_offset + curr_row_size)), ...
            col_ind((col_offset + 1) : (col_offset + curr_col_size)));
        row_offset = row_offset + curr_row_size;
        col_offset = col_offset + curr_col_size;
    end
end

end