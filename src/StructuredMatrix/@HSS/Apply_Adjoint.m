function u = Apply_Adjoint(A, f)
% Apply_Adjoint u = A' * f.

arguments (Input)
    A HSS;
    f (:, :) double;
end

arguments (Output)
    u (:, :) double;
end

if A.leaf_ == 1
    u = A.Amat_' * f;
else
    data = Apply_Adjoint_Upward(A, f);
    u = Apply_Adjoint_Downward(A, data, zeros(A.col_rank_, size(f, 2)));
end

end
