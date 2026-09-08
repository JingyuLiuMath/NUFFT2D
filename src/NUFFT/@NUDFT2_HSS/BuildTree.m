function p = BuildTree(A, x, n_leaf)
% BuildTree Build the binary tree and return its row permutation.

arguments (Input)
    A NUDFT2_HSS;
    x (:, 1) double;
    n_leaf (1, 1) double;
end

arguments (Output)
    p (:, 1) double;
end

A.row_size_ = size(x, 1);
if A.level_ == 0
    A.row_global_size_ = A.row_size_;
end
N = A.col_global_size_;
eta = ceil(A.row_global_size_ / N);
A.max_level_ = A.level_;

if A.col_size_ <= n_leaf || A.row_size_ <= eta * n_leaf
    A.leaf_ = 1;
    p = (1 : A.row_size_)';
else
    A.num_children_ = 2;
    A.children_ = cell(1, A.num_children_);
    first_size = floor(A.col_size_ / 2);
    child_size = [first_size, A.col_size_ - first_size];
    p = zeros(A.row_size_, 1);
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        A.children_{i} = NUDFT2_HSS(N);
        A.children_{i}.col_size_ = child_size(i);
        A.children_{i}.level_ = A.level_ + 1;
        A.children_{i}.row_global_size_ = A.row_global_size_;
        A.children_{i}.row_offset_ = A.row_offset_ + row_offset;
        A.children_{i}.col_offset_ = A.col_offset_ + col_offset;
        I = FindID_ExtendArc(N, A.children_{i}.col_offset_, ...
            A.children_{i}.col_offset_ + child_size(i) - 1, x);
        p_i = A.children_{i}.BuildTree(x(I), n_leaf);
        p(row_offset + (1 : numel(I))) = I(p_i);
        row_offset = row_offset + numel(I);
        col_offset = col_offset + child_size(i);
        A.max_level_ = max(A.max_level_, A.children_{i}.max_level_);
    end
end

end
