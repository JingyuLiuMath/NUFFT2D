function URV_Merge(A)
% URV_Merge

% Jingyu Liu, November 22, 2024.

arguments (Input)
    A HSS;
end

row_size = 0;
col_size = 0;
for i = 1 : A.num_children_
    row_size = row_size + A.children_{i}.URV_row_sk_size_;
    col_size = col_size + A.children_{i}.URV_col_sk_size_;
end

% Merge A.
A.Amat_ = zeros(row_size, col_size);
row_offset = 0;
for i = 1 : A.num_children_
    current_row_size = A.children_{i}.URV_row_sk_size_;
    col_offset = 0;

    for j = 1 : (i - 1)
        current_col_size = A.children_{j}.URV_col_sk_size_;
        A.Amat_((row_offset + 1) : (row_offset + current_row_size), ...
            (col_offset + 1) : (col_offset + current_col_size)) ...
            = A.children_{i}.URV_U_sk_ ...
            * A.Bmat_{i, j} ...
            * A.children_{j}.URV_V_sk_';
        col_offset = col_offset + current_col_size;
    end

    j = i;
    current_col_size = A.children_{j}.URV_col_sk_size_;
    A.Amat_((row_offset + 1) : (row_offset + current_row_size), ...
        (col_offset + 1) : (col_offset + current_col_size)) ...
        = A.children_{i}.URV_A_sk_sk_;
    col_offset = col_offset + current_col_size;

    for j = (i + 1) : A.num_children_
        current_col_size = A.children_{j}.URV_col_sk_size_;
        A.Amat_((row_offset + 1) : (row_offset + current_row_size), ...
            (col_offset + 1) : (col_offset + current_col_size)) ...
            = A.children_{i}.URV_U_sk_ ...
            * A.Bmat_{i, j} ...
            * A.children_{j}.URV_V_sk_';
        col_offset = col_offset + current_col_size;
    end
    row_offset = row_offset + current_row_size;
end

if A.level_ ~= 0
    % Merge U.
    A.Umat_ = zeros(row_size, A.row_rank_);
    row_offset = 0;
    for i = 1 : A.num_children_
        current_row_size = A.children_{i}.URV_row_sk_size_;
        A.Umat_((row_offset + 1) : (row_offset + current_row_size), :) ...
            = A.children_{i}.URV_U_sk_ * A.Rmat_{i};
        row_offset = row_offset + current_row_size;
    end

    % Merge V.
    A.Vmat_ = zeros(col_size, A.col_rank_);
    col_offset = 0;
    for i = 1 : A.num_children_
        current_col_size = A.children_{i}.URV_col_sk_size_;
        A.Vmat_((col_offset + 1) : (col_offset + current_col_size), :) ...
            = A.children_{i}.URV_V_sk_ * A.Wmat_{i};
        col_offset = col_offset + current_col_size;
    end
end

end