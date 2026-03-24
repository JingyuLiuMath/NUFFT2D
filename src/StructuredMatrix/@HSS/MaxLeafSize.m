function [max_row_leaf_size, max_col_leaf_size] = MaxLeafSize(A)
% MaxLeafSize

arguments (Input)
    A HSS;
end

arguments (Output)
    max_row_leaf_size (1, 1) double;
    max_col_leaf_size (1, 1) double;
end

if A.leaf_ == 1
    max_row_leaf_size = A.row_size_;
    max_col_leaf_size = A.col_size_;
else
    max_row_leaf_size = 0;
    max_col_leaf_size = 0;
    for i = 1 : A.num_children_
        [c_max_row_leaf_size, c_max_col_leaf_size] ...
            = A.children_{i}.MaxLeafSize();
        max_row_leaf_size = max(max_row_leaf_size, c_max_row_leaf_size);
        max_col_leaf_size = max(max_col_leaf_size, c_max_col_leaf_size);
    end
end

end