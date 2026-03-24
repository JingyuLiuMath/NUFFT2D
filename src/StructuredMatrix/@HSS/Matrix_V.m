function V = Matrix_V(A)
% Matrix_V

% JingyV Liu, November 20, 2024.

arguments (Input)
    A HSS;
end

arguments (Output)
    V (:, :) double;
end

if A.leaf_ == 1
    V = A.Vmat_;
else
    V = zeros(A.col_size_, A.col_rank_);
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.col_size_;
        V((offset + 1) : (offset + current_size), :) ...
            = A.children_{i}.Matrix_V() * A.Wmat_{i};
        offset = offset + current_size;
    end
end


end