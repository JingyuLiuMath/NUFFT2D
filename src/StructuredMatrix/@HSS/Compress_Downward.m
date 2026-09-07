function Compress_Downward(A, level, tol)
% Top-down SVD stage of HSS recompression, equations (5.9)-(5.16).
% tol is relative to the largest singular value of each local block.

arguments (Input)
    A HSS;
    level (1, 1) double;
    tol (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 0
        P = cell(1, A.num_children_);
        Q = P;
        % Form all sibling bases before changing any coupling matrices.
        for i = 1 : A.num_children_
            Bu = cell(1, A.num_children_);
            Bv = Bu;
            Bu{i} = A.Rmat_{i} * A.Compress_Su_;
            Bv{i} = A.Wmat_{i} * A.Compress_Sv_;
            for j = 1 : A.num_children_
                if i ~= j
                    Bu{j} = A.Bmat_{i, j};
                    Bv{j} = A.Bmat_{j, i}';
                end
            end
            [P{i}, A.children_{i}.Compress_Su_] = MySVDSketch(horzcat(Bu{:}), tol);
            [Q{i}, A.children_{i}.Compress_Sv_] = MySVDSketch(horzcat(Bv{:}), tol);
        end

        for i = 1 : A.num_children_
            for j = 1 : A.num_children_
                if i ~= j
                    A.Bmat_{i, j} = P{i}' * A.Bmat_{i, j} * Q{j};
                end
            end
            A.Rmat_{i} = P{i}' * A.Rmat_{i};
            A.Wmat_{i} = Q{i}' * A.Wmat_{i};
            C = A.children_{i};
            C.row_rank_ = size(P{i}, 2);
            C.col_rank_ = size(Q{i}, 2);
            if C.leaf_ == 1
                C.Umat_ = C.Umat_ * P{i};
                C.Vmat_ = C.Vmat_ * Q{i};
            else
                % Pass the parent basis changes down, as in (5.10).
                for j = 1 : C.num_children_
                    C.Rmat_{j} = C.Rmat_{j} * P{i};
                    C.Wmat_{j} = C.Wmat_{j} * Q{i};
                end
            end
        end
    end

    % Clear, including leaf nodes.
    A.Compress_Su_ = [];
    A.Compress_Sv_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.Compress_Downward(level, tol);
    end
end

end
