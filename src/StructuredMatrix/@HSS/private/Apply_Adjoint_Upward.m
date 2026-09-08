function data = Apply_Adjoint_Upward(A, f)
% Project onto the nested U bases and apply the adjoint leaf blocks.

data.children = cell(1, A.num_children_);
if A.leaf_ == 1
    data.fhvec = A.Umat_' * f;
    data.uvec = A.Amat_' * f;
else
    data.fhvec = zeros(A.row_rank_, size(f, 2));
    row_offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.row_size_;
        data.children{i} = Apply_Adjoint_Upward(A.children_{i}, ...
            f((row_offset + 1) : (row_offset + current_size), :));
        if A.level_ ~= 0
            data.fhvec = data.fhvec ...
                + A.children_{i}.Rmat_' * data.children{i}.fhvec;
        end
        row_offset = row_offset + current_size;
    end
end

end
