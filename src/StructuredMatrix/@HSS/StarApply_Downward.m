function StarApply_Downward(A, level)
% StarApply_Downward

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.uvec_ = A.Amat_' * A.fvec_ +  A.Vmat_ * A.uhvec_;
    else
        for i = 1 : A.num_children_
            A.children_{i}.uhvec_ = A.Wmat_{i} * A.uhvec_;
        end

        for i = 1 : A.num_children_
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                A.children_{i}.uhvec_ = A.children_{i}.uhvec_ ...
                    + A.Bmat_{j, i}' * A.children_{j}.fhvec_;
            end
        end
    end

    % Clear.
    A.fvec_ = [];
    A.fhvec_ = [];
    A.uhvec_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.StarApply_Downward(level);
    end
end

end