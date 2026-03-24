function StarApply_Upward(A, level)
% StarApply_Upward

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.fhvec_ = A.Umat_' * A.fvec_;
    else
        A.fhvec_ = zeros(A.row_rank_, A.vec_col_size_);
        for i = 1 : A.num_children_
            A.fhvec_ = A.fhvec_ + A.Rmat_{i}' * A.children_{i}.fhvec_;
        end
    end
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.StarApply_Upward(level);
    end
end

end