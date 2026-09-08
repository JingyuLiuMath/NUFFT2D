function data = URV_Merge(A, children, is_root)
% Assemble the reduced diagonal block and nested bases from child workspaces.

row_size = 0;
col_size = 0;
for i = 1 : A.num_children_
    row_size = row_size + size(children{i}.Amat, 1);
    col_size = col_size + size(children{i}.Amat, 2);
end

data.Amat = zeros(row_size, col_size);
row_offset = 0;
for i = 1 : A.num_children_
    current_row_size = size(children{i}.Amat, 1);
    col_offset = 0;
    for j = 1 : A.num_children_
        current_col_size = size(children{j}.Amat, 2);
        if i == j
            data.Amat((row_offset + 1) : (row_offset + current_row_size), ...
                (col_offset + 1) : (col_offset + current_col_size)) ...
                = children{i}.Amat;
        else
            data.Amat((row_offset + 1) : (row_offset + current_row_size), ...
                (col_offset + 1) : (col_offset + current_col_size)) ...
                = children{i}.Umat * A.Bmat_{i, j} * children{j}.Vmat';
        end
        col_offset = col_offset + current_col_size;
    end
    row_offset = row_offset + current_row_size;
end

if is_root
    return;
end

data.Umat = zeros(row_size, A.row_rank_);
data.Vmat = zeros(col_size, A.col_rank_);
row_offset = 0;
col_offset = 0;
for i = 1 : A.num_children_
    current_row_size = size(children{i}.Amat, 1);
    current_col_size = size(children{i}.Amat, 2);
    data.Umat((row_offset + 1) : (row_offset + current_row_size), :) ...
        = children{i}.Umat * A.children_{i}.Rmat_;
    data.Vmat((col_offset + 1) : (col_offset + current_col_size), :) ...
        = children{i}.Vmat * A.children_{i}.Wmat_;
    row_offset = row_offset + current_row_size;
    col_offset = col_offset + current_col_size;
end

end
