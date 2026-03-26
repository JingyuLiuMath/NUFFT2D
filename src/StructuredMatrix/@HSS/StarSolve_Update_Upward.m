function StarSolve_Update_Upward(A, level, max_level)
% StarSolve_Update_Upward

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
    max_level (1, 1) double;
end

if A.level_ == level
    if A.level_ == max_level
        A.fhvec_ = A.URV_U_re_' * A.URV_f_re_;
    else
        A.fhvec_ = zeros(A.row_rank_, A.vec_col_size_);
        for i = 1 : A.num_children_
            A.fhvec_ = A.fhvec_ + A.Rmat_{i}' * A.children_{i}.fhvec_;
        end
    end
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.StarSolve_Update_Upward(level, max_level);
    end
end

end