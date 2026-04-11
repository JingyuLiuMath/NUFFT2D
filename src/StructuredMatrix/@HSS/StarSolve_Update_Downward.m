function StarSolve_Update_Downward(A, level, max_level)
% StarSolve_Update_Downward

% Jingyu Liu, December 8, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
    max_level (1, 1) double;
end

if A.level_ == level
    if A.level_ == max_level || A.urv_leaf_ == 1
        if isempty(A.URV_u_sk_)
            A.URV_u_sk_ = -A.URV_V_sk_ * A.uhvec_;
        else
            A.URV_u_sk_ = A.URV_u_sk_ ...
                - A.URV_A_re_sk_' * A.URV_f_re_ -  A.URV_V_sk_ * A.uhvec_;
        end
    elseif A.leaf_ == 0
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
    A.fhvec_ = [];
    A.uhvec_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.StarSolve_Update_Downward(level, max_level);
    end
end

end