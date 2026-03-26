function URV_StarSolve_SplitVector(A)
% URV_StarSolve_SplitVector

% Jingyu Liu, December 9, 2024.

arguments (Input)
    A HSS;
end

A.urv_leaf_ = A.leaf_;

offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.URV_row_sk_size_;
    A.children_{i}.URV_f_sk_ ...
        = A.fvec_((offset + 1) : (offset + current_size), :);
    offset = offset + current_size;
end

% Clear.
A.fvec_ = [];

end