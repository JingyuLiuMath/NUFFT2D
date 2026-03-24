function URV_Solve_SplitVector(A)
% URV_Solve_SplitVector

% Jingyu Liu, November 23, 2024.

arguments (Input)
    A HSS;
end

offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.URV_col_sk_size_;
    A.children_{i}.URV_u_sk_ ...
        = A.uvec_((offset + 1) : (offset + current_size), :);
    offset = offset + current_size;
end

if A.level_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.fhvec_ = 0;
    end
else
    for i = 1 : A.num_children_
        A.children_{i}.fhvec_ = A.Rmat_{i} * A.fhvec_;
    end
end

% Compute fh;
for i = 1 : A.num_children_
    for j = [1 : (i - 1), (i + 1) : A.num_children_]
        A.children_{i}.fhvec_ ...
            = A.children_{i}.fhvec_ ...
            + A.Bmat_{i, j} * ...
            (A.children_{j}.URV_V_sk_' * A.children_{j}.URV_u_sk_);
    end
end

% Clear.
A.uvec_ = [];
A.fhvec_ = [];

end