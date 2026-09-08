function V = Matrix_V(A)
% Matrix_V Expand the nested V basis.

arguments (Input)
    A HSS;
end

arguments (Output)
    V (:, :) double;
end

if A.level_ == 0
    V = zeros(A.col_size_, 0);
elseif A.leaf_ == 1
    V = A.Vmat_;
else
    V = zeros(A.col_size_, A.col_rank_);
    col_offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.col_size_;
        V((col_offset + 1) : (col_offset + current_size), :) ...
            = A.children_{i}.Matrix_V() * A.children_{i}.Wmat_;
        col_offset = col_offset + current_size;
    end
end

end
