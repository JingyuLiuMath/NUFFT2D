function StarSolve_Update_Root(A)
% StarSolve_Update_Root

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A HSS;
end

if A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.uhvec_ ...
            = zeros(A.children_{i}.col_rank_, A.vec_col_size_);
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            A.children_{i}.uhvec_ = A.children_{i}.uhvec_ ...
                + A.Bmat_{j, i}' * A.children_{j}.fhvec_;
        end
    end
end


end