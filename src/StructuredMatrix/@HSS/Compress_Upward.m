function Compress_Upward(A, level)
% Bottom-up QR stage of HSS recompression, equations (5.5)-(5.8).

arguments (Input)
    A HSS;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        U = A.Umat_;
        V = A.Vmat_;
    else
        row_size = 0;
        col_size = 0;
        for i = 1 : A.num_children_
            row_size = row_size + A.children_{i}.row_rank_;
            col_size = col_size + A.children_{i}.col_rank_;
        end
        U = zeros(row_size, A.row_rank_);
        V = zeros(col_size, A.col_rank_);
        row_offset = 0;
        col_offset = 0;
        for i = 1 : A.num_children_
            U(row_offset + (1 : A.children_{i}.row_rank_), :) ...
                = A.children_{i}.Compress_Tu_ * A.Rmat_{i};
            V(col_offset + (1 : A.children_{i}.col_rank_), :) ...
                = A.children_{i}.Compress_Tv_ * A.Wmat_{i};
            row_offset = row_offset + A.children_{i}.row_rank_;
            col_offset = col_offset + A.children_{i}.col_rank_;
        end
        for i = 1 : A.num_children_
            for j = 1 : A.num_children_
                if i ~= j
                    A.Bmat_{i, j} = A.children_{i}.Compress_Tu_ ...
                        * A.Bmat_{i, j} * A.children_{j}.Compress_Tv_';
                else
                    A.Bmat_{i, j} = [];
                end
            end
        end
        % Clear after all sibling couplings have used the child transforms.
        for i = 1 : A.num_children_
            A.children_{i}.Compress_Tu_ = [];
            A.children_{i}.Compress_Tv_ = [];
        end
    end

    [U, A.Compress_Tu_] = qr(U, "econ");
    [V, A.Compress_Tv_] = qr(V, "econ");
    A.row_rank_ = size(U, 2);
    A.col_rank_ = size(V, 2);
    if A.leaf_ == 1
        A.Umat_ = U;
        A.Vmat_ = V;
    else
        row_offset = 0;
        col_offset = 0;
        for i = 1 : A.num_children_
            A.Rmat_{i} = U(row_offset + (1 : A.children_{i}.row_rank_), :);
            A.Wmat_{i} = V(col_offset + (1 : A.children_{i}.col_rank_), :);
            row_offset = row_offset + A.children_{i}.row_rank_;
            col_offset = col_offset + A.children_{i}.col_rank_;
        end
    end
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.Compress_Upward(level);
    end
end

end
