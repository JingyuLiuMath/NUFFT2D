function u = Apply_Adjoint_Downward(A, data, uhvec)
% Propagate incoming coefficients and add adjoint sibling interactions.

if A.leaf_ == 1
    u = data.uvec + A.Vmat_ * uhvec;
else
    u = zeros(A.col_size_, size(uhvec, 2));
    col_offset = 0;
    for i = 1 : A.num_children_
        if A.level_ == 0
            uhvec_i = zeros(A.children_{i}.col_rank_, size(uhvec, 2));
        else
            uhvec_i = A.children_{i}.Wmat_ * uhvec;
        end
        for j = 1 : A.num_children_
            if i ~= j
                uhvec_i = uhvec_i + A.Bmat_{j, i}' * data.children{j}.fhvec;
            end
        end
        current_size = A.children_{i}.col_size_;
        u((col_offset + 1) : (col_offset + current_size), :) ...
            = Apply_Adjoint_Downward(A.children_{i}, data.children{i}, uhvec_i);
        col_offset = col_offset + current_size;
    end
end

end
