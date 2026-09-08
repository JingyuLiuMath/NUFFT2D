function Recompress_Downward(A, Su, Sv, tol, is_root)
% Compress sibling interactions together with the inherited weights.

if A.leaf_ == 1
    return;
end

row_size = 0;
col_size = 0;
for i = 1 : A.num_children_
    row_size = row_size + A.children_{i}.row_rank_;
    col_size = col_size + A.children_{i}.col_rank_;
end

data.children = cell(1, A.num_children_);
% Form every sibling basis before changing the coupling matrices.
for i = 1 : A.num_children_
    Bu = zeros(A.children_{i}.row_rank_, ...
        col_size - A.children_{i}.col_rank_ + size(Su, 2));
    Bv = zeros(A.children_{i}.col_rank_, ...
        row_size - A.children_{i}.row_rank_ + size(Sv, 2));
    row_offset = 0;
    col_offset = 0;
    for j = 1 : A.num_children_
        if i == j
            if is_root
                continue;
            end
            current_row_size = size(Su, 2);
            current_col_size = size(Sv, 2);
            Bu(:, row_offset + (1 : current_row_size)) = A.children_{i}.Rmat_ * Su;
            Bv(:, col_offset + (1 : current_col_size)) = A.children_{i}.Wmat_ * Sv;
        else
            current_row_size = A.children_{j}.col_rank_;
            current_col_size = A.children_{j}.row_rank_;
            Bu(:, row_offset + (1 : current_row_size)) = A.Bmat_{i, j};
            Bv(:, col_offset + (1 : current_col_size)) = A.Bmat_{j, i}';
        end
        row_offset = row_offset + current_row_size;
        col_offset = col_offset + current_col_size;
    end
    [data.children{i}.P, data.children{i}.Su] = MySVDSketch(Bu, tol);
    [data.children{i}.Q, data.children{i}.Sv] = MySVDSketch(Bv, tol);
end

for i = 1 : A.num_children_
    for j = 1 : A.num_children_
        if i ~= j
            A.Bmat_{i, j} = data.children{i}.P' * A.Bmat_{i, j} * data.children{j}.Q;
        end
    end
    if ~is_root
        A.children_{i}.Rmat_ = data.children{i}.P' * A.children_{i}.Rmat_;
        A.children_{i}.Wmat_ = data.children{i}.Q' * A.children_{i}.Wmat_;
    end
    A.children_{i}.row_rank_ = size(data.children{i}.P, 2);
    A.children_{i}.col_rank_ = size(data.children{i}.Q, 2);
    if A.children_{i}.leaf_ == 1
        A.children_{i}.Umat_ = A.children_{i}.Umat_ * data.children{i}.P;
        A.children_{i}.Vmat_ = A.children_{i}.Vmat_ * data.children{i}.Q;
    else
        for j = 1 : A.children_{i}.num_children_
            A.children_{i}.children_{j}.Rmat_ ...
                = A.children_{i}.children_{j}.Rmat_ * data.children{i}.P;
            A.children_{i}.children_{j}.Wmat_ ...
                = A.children_{i}.children_{j}.Wmat_ * data.children{i}.Q;
        end
    end
end

for i = 1 : A.num_children_
    Recompress_Downward(A.children_{i}, data.children{i}.Su, ...
        data.children{i}.Sv, tol, false);
end

end
