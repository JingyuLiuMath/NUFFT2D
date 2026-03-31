function URV_Solve_Root(A)
% URV_Solve_Root

% Jingyu Liu, November 23, 2024.

arguments (Input)
    A HSS;
end

if A.leaf_ == 0
    A.URV_Solve_MergeVector();
end

m = size(A.URV_A_re_re_, 1);
n = size(A.URV_A_re_re_, 2);

A.fvec_ = A.URV_P_' * A.fvec_;
if m >= n
    A.uvec_ = A.URV_A_re_re_(1 : n, :) \ A.fvec_(1 : n, :);
else
    A.uvec = [A.URV_A_re_re_(:, 1 : m) \ A.fvec_; ...
        zeros(n - m, A.vec_col_size_)];
end

% Clear.
A.fvec_ = [];

if A.leaf_ == 0
    A.URV_Solve_SplitVector();
end

end