function U = Matrix_U(A)
% Matrix_U Expand the nested U basis.

arguments (Input)
    A HSS;
end

arguments (Output)
    U (:, :) double;
end

if A.level_ == 0
    U = zeros(A.row_size_, 0);
elseif A.leaf_ == 1
    U = A.Umat_;
else
    U = zeros(A.row_size_, A.row_rank_);
    row_offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.row_size_;
        U((row_offset + 1) : (row_offset + current_size), :) ...
            = A.children_{i}.Matrix_U() * A.children_{i}.Rmat_;
        row_offset = row_offset + current_size;
    end
end

end
