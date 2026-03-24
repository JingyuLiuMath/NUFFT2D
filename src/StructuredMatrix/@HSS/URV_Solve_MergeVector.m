function URV_Solve_MergeVector(A)
% URV_Solve_MergeVector

% Jingyu Liu, November 23, 2024.

arguments (Input)
    A HSS;
end

row_size = 0;
for i = 1 : A.num_children_
    row_size = row_size + A.children_{i}.URV_row_sk_size_;
end

A.fvec_ = zeros(row_size, A.vec_col_size_);
offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.URV_row_sk_size_;
    A.fvec_((offset + 1) : (offset + current_size), :) ...
        = A.children_{i}.URV_f_sk_;
    offset = offset + current_size;

    % Clear.
    A.children_{i}.URV_f_sk_ = [];
end

end