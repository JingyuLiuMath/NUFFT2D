function [row_leaf, col_leaf, num_leaves] ...
    = LeafRows(A, row_leaf, col_leaf, num_leaves)

if nargin == 1
    row_leaf = zeros(A.row_size_, 1);
    col_leaf = zeros(A.col_size_, 1);
    num_leaves = 0;
end

if A.leaf_ == 1
    num_leaves = num_leaves + 1;
    row_leaf(A.row_offset_ + (1 : A.row_size_)) = num_leaves;
    col_leaf(A.col_offset_ + (1 : A.col_size_)) = num_leaves;
else
    for i = 1 : A.num_children_
        [row_leaf, col_leaf, num_leaves] = LeafRows(...
            A.children_{i}, row_leaf, col_leaf, num_leaves);
    end
end

end
