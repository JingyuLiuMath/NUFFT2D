function [p, q, row_index] = BuildTree_FaceSplitting(A, Ax, Ay, px, py)
% Build the 2D tree and return row indices in the two 1D trees.

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    px (:, 1) double;
    py (:, 1) double;
end

M = size(px, 1);
x_inv_perm = zeros(M, 1);
y_inv_perm = zeros(M, 1);
x_inv_perm(px) = (1 : M)';
y_inv_perm(py) = (1 : M)';

% Count the rows in each pair of 1D leaves.
[row_leaf_x, data.col_leaf_x, data.num_leaves_x] = LeafRows(Ax);
[row_leaf_y, data.col_leaf_y, data.num_leaves_y] = LeafRows(Ay);
leaf_pair = row_leaf_x(x_inv_perm) ...
    + data.num_leaves_x * (row_leaf_y(y_inv_perm) - 1);
num_leaves = data.num_leaves_x * data.num_leaves_y;
data.row_counts = accumarray(leaf_pair, ones(M, 1), [num_leaves, 1]);

% Tree and leaf offsets.
leaf_row_offset = zeros(num_leaves, 1);
p = zeros(M, 1);
q = zeros(Ax.col_size_ * Ay.col_size_, 1);
A.row_global_size_ = M;
[q, leaf_row_offset] = BuildNode(A, Ax, Ay, data, q, leaf_row_offset);

% Preserve the original row order within each 2D leaf.
for i = 1 : M
    k = leaf_pair(i);
    leaf_row_offset(k) = leaf_row_offset(k) + 1;
    p(leaf_row_offset(k)) = i;
end
row_index.x = x_inv_perm(p);
row_index.y = y_inv_perm(p);

end

function [row_leaf, col_leaf, num_leaves] = LeafRows(A)
% Leaf numbers for source rows and leaf column starts.

row_leaf = zeros(A.row_size_, 1);
col_leaf = zeros(A.col_size_, 1);
[row_leaf, col_leaf, num_leaves] = LabelLeaves(A, row_leaf, col_leaf, 0);

end

function [row_leaf, col_leaf, num_leaves] = LabelLeaves(A, row_leaf, col_leaf, num_leaves)

if A.leaf_ == 1
    num_leaves = num_leaves + 1;
    row_leaf(A.row_offset_ + (1 : A.row_size_)) = num_leaves;
    col_leaf(A.col_offset_ + 1) = num_leaves;
else
    for i = 1 : A.num_children_
        [row_leaf, col_leaf, num_leaves] ...
            = LabelLeaves(A.children_{i}, row_leaf, col_leaf, num_leaves);
    end
end

end

function [q, leaf_row_offset] = BuildNode(A, Ax, Ay, data, q, leaf_row_offset)

A.col_size_ = Ax.col_size_ * Ay.col_size_;
A.max_level_ = A.level_;
A.leaf_ = double(Ax.leaf_ == 1 && Ay.leaf_ == 1);
if A.leaf_ == 1
    k = data.col_leaf_x(Ax.col_offset_ + 1) ...
        + data.num_leaves_x * (data.col_leaf_y(Ay.col_offset_ + 1) - 1);
    A.row_size_ = data.row_counts(k);
    leaf_row_offset(k) = A.row_offset_;
    col_x = (Ax.col_offset_ : (Ax.col_offset_ + Ax.col_size_ - 1))';
    col_y = Ay.col_offset_ : (Ay.col_offset_ + Ay.col_size_ - 1);
    q(A.col_offset_ + (1 : A.col_size_)) ...
        = reshape(col_x + Ax.col_global_size_ * col_y + 1, [], 1);
    return;
end

num_children_x = max(Ax.num_children_, 1);
num_children_y = max(Ay.num_children_, 1);
A.num_children_ = num_children_x * num_children_y;
A.children_ = cell(1, A.num_children_);
row_offset = 0;
col_offset = 0;
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    A.children_{i} = NUDFT2_2D_HSS(Ax.col_global_size_, Ay.col_global_size_);
    A.children_{i}.level_ = A.level_ + 1;
    A.children_{i}.row_global_size_ = A.row_global_size_;
    A.children_{i}.row_offset_ = A.row_offset_ + row_offset;
    A.children_{i}.col_offset_ = A.col_offset_ + col_offset;
    if Ax.leaf_ == 1
        [q, leaf_row_offset] = BuildNode(A.children_{i}, Ax, Ay.children_{iy}, ...
            data, q, leaf_row_offset);
    elseif Ay.leaf_ == 1
        [q, leaf_row_offset] = BuildNode(A.children_{i}, Ax.children_{ix}, Ay, ...
            data, q, leaf_row_offset);
    else
        [q, leaf_row_offset] = BuildNode(A.children_{i}, Ax.children_{ix}, Ay.children_{iy}, ...
            data, q, leaf_row_offset);
    end
    row_offset = row_offset + A.children_{i}.row_size_;
    col_offset = col_offset + A.children_{i}.col_size_;
    A.max_level_ = max(A.max_level_, A.children_{i}.max_level_);
end
A.row_size_ = row_offset;

end
