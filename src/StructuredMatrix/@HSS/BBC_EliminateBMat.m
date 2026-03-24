function BBC_EliminateBMat(A, level)
% BBC_EliminateBMat

% Jingyu Liu, December 6, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.Amat_ = A.Amat_ + A.Umat_ * A.BBC_A_ * A.Vmat_';
        A.BBC_A_ = [];
    else
        % Update B.
        for i = 1 : A.num_children_
            for j = 1 : A.num_children_
                A.Bmat_{i, j} = A.Bmat_{i, j} ...
                    + A.Rmat_{i} * A.BBC_A_ * A.Wmat_{j}';
            end
        end
        A.BBC_A_ = [];

        % Splist B.
        for i = 1 : A.num_children_
            A.children_{i}.BBC_A_ = A.Bmat_{i, i};
            A.Bmat_{i, i} = [];
        end
    end
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.BBC_EliminateBMat(level);
    end
end

end