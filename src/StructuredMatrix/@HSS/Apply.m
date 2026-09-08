function f = Apply(A, u)
% Apply f = A * u.

arguments (Input)
    A HSS;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

if A.leaf_ == 1
    f = A.Amat_ * u;
else
    data = Apply_Upward(A, u);
    f = Apply_Downward(A, data, zeros(A.row_rank_, size(u, 2)));
end

end
