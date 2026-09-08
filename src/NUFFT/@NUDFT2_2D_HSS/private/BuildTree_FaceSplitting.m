function [p, q, index] = BuildTree_FaceSplitting(A, Ax, Ay, px, py)
% Build the product tree and group original rows once by their leaf pair.

M = size(px, 1);
index.posx = zeros(M, 1);
index.posy = zeros(M, 1);
index.posx(px) = (1 : M)';
index.posy(py) = (1 : M)';
[leaf_x, data.start_x, data.num_x] = LeafRows(Ax);
[leaf_y, data.start_y, num_y] = LeafRows(Ay);
bucket = leaf_x(index.posx) + data.num_x * (leaf_y(index.posy) - 1);
data.counts = accumarray(bucket, ones(M, 1), [data.num_x * num_y, 1]);
bucket_offset = zeros(data.num_x * num_y, 1);
p = zeros(M, 1);
q = zeros(Ax.col_size_ * Ay.col_size_, 1);
A.row_global_size_ = M;
[q, bucket_offset] = BuildNode(A, Ax, Ay, data, q, bucket_offset);

% Original row order makes every bucket agree with intersect(..., "sorted").
for sample = 1 : M
    bucket_offset(bucket(sample)) = bucket_offset(bucket(sample)) + 1;
    p(bucket_offset(bucket(sample))) = sample;
end

end

function [q, bucket_offset] = BuildNode(A, Ax, Ay, data, q, bucket_offset)

A.col_size_ = Ax.col_size_ * Ay.col_size_;
A.max_level_ = A.level_;
A.leaf_ = double(Ax.leaf_ == 1 && Ay.leaf_ == 1);
if A.leaf_ == 1
    k = data.start_x(Ax.col_offset_ + 1) ...
        + data.num_x * (data.start_y(Ay.col_offset_ + 1) - 1);
    A.row_size_ = data.counts(k);
    bucket_offset(k) = A.row_offset_;
    q(A.col_offset_ + (1 : A.col_size_)) = reshape( ...
        (Ax.col_offset_ : (Ax.col_offset_ + Ax.col_size_ - 1))' ...
        + Ax.col_global_size_ ...
        * (Ay.col_offset_ : (Ay.col_offset_ + Ay.col_size_ - 1)) + 1, [], 1);
    return;
end

child_x = max(Ax.num_children_, 1);
child_y = max(Ay.num_children_, 1);
A.num_children_ = child_x * child_y;
A.children_ = cell(1, A.num_children_);
row_offset = 0;
col_offset = 0;
for i = 1 : A.num_children_
    ix = mod(i - 1, child_x) + 1;
    iy = floor((i - 1) / child_x) + 1;
    A.children_{i} = NUDFT2_2D_HSS(Ax.col_global_size_, Ay.col_global_size_);
    A.children_{i}.level_ = A.level_ + 1;
    A.children_{i}.row_global_size_ = A.row_global_size_;
    A.children_{i}.row_offset_ = A.row_offset_ + row_offset;
    A.children_{i}.col_offset_ = A.col_offset_ + col_offset;
    if Ax.leaf_ == 1
        [q, bucket_offset] = BuildNode(A.children_{i}, Ax, Ay.children_{iy}, ...
            data, q, bucket_offset);
    elseif Ay.leaf_ == 1
        [q, bucket_offset] = BuildNode(A.children_{i}, Ax.children_{ix}, Ay, ...
            data, q, bucket_offset);
    else
        [q, bucket_offset] = BuildNode(A.children_{i}, Ax.children_{ix}, Ay.children_{iy}, ...
            data, q, bucket_offset);
    end
    row_offset = row_offset + A.children_{i}.row_size_;
    col_offset = col_offset + A.children_{i}.col_size_;
    A.max_level_ = max(A.max_level_, A.children_{i}.max_level_);
end
A.row_size_ = row_offset;

end

function [leaf, start, num_leaves] = LeafRows(A)
% Label source rows and column starts without copying indices up the tree.

leaf = zeros(A.row_size_, 1);
start = zeros(A.col_size_, 1);
[leaf, start, num_leaves] = Visit(A, leaf, start, 0);

end

function [leaf, start, num_leaves] = Visit(A, leaf, start, num_leaves)

if A.leaf_ == 1
    num_leaves = num_leaves + 1;
    leaf(A.row_offset_ + (1 : A.row_size_)) = num_leaves;
    start(A.col_offset_ + 1) = num_leaves;
else
    for i = 1 : A.num_children_
        [leaf, start, num_leaves] = Visit(A.children_{i}, leaf, start, num_leaves);
    end
end

end
