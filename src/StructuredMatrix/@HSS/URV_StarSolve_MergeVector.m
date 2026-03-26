function URV_StarSolve_MergeVector(A)
% URV_StarSolve_MergeVector

% Jingyu Liu, December 9, 2024.

arguments (Input)
    A HSS;
end

A.urv_leaf_ = 1;

col_size = 0;
for i = 1 : A.num_children_
    col_size = col_size + A.children_{i}.URV_col_sk_size_;
end

A.uvec_ = zeros(col_size, A.vec_col_size_);
offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.URV_col_sk_size_;
    A.uvec_((offset + 1) : (offset + current_size), :) ...
        = A.children_{i}.URV_u_sk_;
    offset = offset + current_size;

    % Clear.
    A.children_{i}.URV_u_sk_ = [];
end
end