function Apply_Root(A)
% Apply_Root

% Jingyu Liu, November 18, 2024.

arguments (Input)
    A HSS;
end

if A.leaf_ == 1
    A.fvec_ = A.Amat_ * A.uvec_;
else
    for i = 1 : A.num_children_
        A.children_{i}.fhvec_ ...
            = zeros(A.children_{i}.row_rank_, A.vec_col_size_);
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            A.children_{i}.fhvec_ = A.children_{i}.fhvec_ ...
                + A.Bmat_{i, j} * A.children_{j}.uhvec_;
        end
    end
end


end