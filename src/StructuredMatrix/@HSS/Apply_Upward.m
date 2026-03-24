function Apply_Upward(A, level)
% Apply_Upward

% Jingyu Liu, November 18, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.uhvec_ = A.Vmat_' * A.uvec_;
    else
        A.uhvec_ = zeros(A.col_rank_, A.vec_col_size_);
        for i = 1 : A.num_children_
            A.uhvec_ = A.uhvec_ + A.Wmat_{i}' * A.children_{i}.uhvec_;
        end
    end
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.Apply_Upward(level);
    end
end

end