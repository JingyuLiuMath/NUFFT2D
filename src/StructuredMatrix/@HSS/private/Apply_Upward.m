function data = Apply_Upward(A, u)
% Project onto the nested V bases and apply the dense leaf blocks.

data.children = cell(1, A.num_children_);
if A.leaf_ == 1
    data.uhvec = A.Vmat_' * u;
    data.fvec = A.Amat_ * u;
else
    data.uhvec = zeros(A.col_rank_, size(u, 2));
    col_offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.col_size_;
        data.children{i} = Apply_Upward(A.children_{i}, ...
            u((col_offset + 1) : (col_offset + current_size), :));
        if A.level_ ~= 0
            data.uhvec = data.uhvec ...
                + A.children_{i}.Wmat_' * data.children{i}.uhvec;
        end
        col_offset = col_offset + current_size;
    end
end

end
