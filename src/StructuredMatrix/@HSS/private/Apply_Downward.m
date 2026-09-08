function f = Apply_Downward(A, data, fhvec)
% Propagate incoming coefficients and add sibling interactions.

if A.leaf_ == 1
    f = data.fvec + A.Umat_ * fhvec;
else
    f = zeros(A.row_size_, size(fhvec, 2));
    row_offset = 0;
    for i = 1 : A.num_children_
        if A.level_ == 0
            fhvec_i = zeros(A.children_{i}.row_rank_, size(fhvec, 2));
        else
            fhvec_i = A.children_{i}.Rmat_ * fhvec;
        end
        for j = 1 : A.num_children_
            if i ~= j
                fhvec_i = fhvec_i + A.Bmat_{i, j} * data.children{j}.uhvec;
            end
        end
        current_size = A.children_{i}.row_size_;
        f((row_offset + 1) : (row_offset + current_size), :) ...
            = Apply_Downward(A.children_{i}, data.children{i}, fhvec_i);
        row_offset = row_offset + current_size;
    end
end

end
