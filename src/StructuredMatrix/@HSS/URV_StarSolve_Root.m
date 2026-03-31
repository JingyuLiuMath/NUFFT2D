function URV_StarSolve_Root(A)
% URV_StarSolve_Root

% Jingyu Liu, December 9, 2024.

arguments (Input)
    A HSS;
end

if A.leaf_ == 0
    A.URV_StarSolve_MergeVector();
end

m = size(A.URV_A_re_re_, 1);
n = size(A.URV_A_re_re_, 2);

if m >= n
    A.fvec_ = [A.URV_A_re_re_(1 : n, :)' \ A.uvec_; ...
        zeros(m - n, A.vec_col_size_)];
else
    A.fvec_ = A.URV_A_re_re_(:, 1 : m)' \ A.uvec_(1 : m, :);
end
A.fvec_ = A.URV_P_ * A.fvec_;

% Clear.
A.uvec_ = [];

if A.leaf_ == 0
    A.URV_StarSolve_SplitVector();
end

end