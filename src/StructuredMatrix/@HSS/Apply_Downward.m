function Apply_Downward(A, level)
% Apply_Downward

% Jingyu Liu, November 18, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.fvec_ = A.Amat_ * A.uvec_ +  A.Umat_ * A.fhvec_;
    else
        for i = 1 : A.num_children_
            A.children_{i}.fhvec_ = A.Rmat_{i} * A.fhvec_;
        end

        for i = 1 : A.num_children_
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                A.children_{i}.fhvec_ = A.children_{i}.fhvec_ ...
                    + A.Bmat_{i, j} * A.children_{j}.uhvec_;
            end
        end
    end

    % Clear.
    A.uvec_ = [];
    A.uhvec_ = [];
    A.fhvec_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.Apply_Downward(level);
    end
end

end