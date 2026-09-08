function data = Recompress_Upward(A, is_root)
% Orthogonalize nested bases and pass the QR factors to the parent.

data = struct;
if A.leaf_ == 1
    U = A.Umat_;
    V = A.Vmat_;
else
    children_data = cell(1, A.num_children_);
    row_size = 0;
    col_size = 0;
    for i = 1 : A.num_children_
        children_data{i} = Recompress_Upward(A.children_{i}, false);
        row_size = row_size + A.children_{i}.row_rank_;
        col_size = col_size + A.children_{i}.col_rank_;
    end

    if ~is_root
        U = zeros(row_size, A.row_rank_);
        V = zeros(col_size, A.col_rank_);
        row_offset = 0;
        col_offset = 0;
        for i = 1 : A.num_children_
            U(row_offset + (1 : A.children_{i}.row_rank_), :) ...
                = children_data{i}.Tu * A.children_{i}.Rmat_;
            V(col_offset + (1 : A.children_{i}.col_rank_), :) ...
                = children_data{i}.Tv * A.children_{i}.Wmat_;
            row_offset = row_offset + A.children_{i}.row_rank_;
            col_offset = col_offset + A.children_{i}.col_rank_;
        end
    end

    for i = 1 : A.num_children_
        for j = 1 : A.num_children_
            if i ~= j
                A.Bmat_{i, j} = children_data{i}.Tu ...
                    * A.Bmat_{i, j} * children_data{j}.Tv';
            else
                A.Bmat_{i, j} = [];
            end
        end
    end
end

if is_root
    return;
end

[U, data.Tu] = qr(U, "econ");
[V, data.Tv] = qr(V, "econ");
A.row_rank_ = size(U, 2);
A.col_rank_ = size(V, 2);
if A.leaf_ == 1
    A.Umat_ = U;
    A.Vmat_ = V;
else
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        A.children_{i}.Rmat_ = U(row_offset + (1 : A.children_{i}.row_rank_), :);
        A.children_{i}.Wmat_ = V(col_offset + (1 : A.children_{i}.col_rank_), :);
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end
end

end
