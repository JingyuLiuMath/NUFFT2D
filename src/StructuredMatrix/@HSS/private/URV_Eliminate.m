function [Fac, data] = URV_Eliminate(A, is_root)
% Recursively eliminate nodes and pass reduced generators to the parent.
% Each non-root node requires m >= n - s, as in the original URV algorithm.

Fac.leaf = A.leaf_;
Fac.num_children = A.num_children_;
Fac.row_size = A.row_size_;
Fac.col_size = A.col_size_;
Fac.row_rank = A.row_rank_;
Fac.col_rank = A.col_rank_;
Fac.children = cell(1, A.num_children_);
Fac.Bmat = A.Bmat_;
Fac.Rmat = A.Rmat_;
Fac.Wmat = A.Wmat_;
Fac.Omega = [];
if is_root
    Fac.Rmat = [];
    Fac.Wmat = [];
end

if A.leaf_ == 1
    data.Amat = A.Amat_;
    data.Umat = A.Umat_;
    data.Vmat = A.Vmat_;
else
    data.children = cell(1, A.num_children_);
    for i = 1 : A.num_children_
        [Fac.children{i}, data.children{i}] ...
            = URV_Eliminate(A.children_{i}, false);
    end
    data = URV_Merge(A, data.children, is_root);
end

if is_root
    return;
end

m = size(data.Amat, 1);
n = size(data.Amat, 2);
k = A.row_rank_;
s = A.col_rank_;

% Reduce rows using the joint range of U and the diagonal block.
if m >= k + n
    Z = zeros(m, k + n);
    Z(:, 1 : k) = data.Umat;
    Z(:, (k + 1) : (k + n)) = data.Amat;
    [Fac.Omega, Z] = qr(Z, "econ");
    m = k + n;
    data.Umat = Z(1 : m, 1 : k);
    data.Amat = Z(1 : m, (k + 1) : (k + n));
end

% Zero the redundant rows of the transformed V basis.
[Fac.Q, R] = qr(data.Vmat);
reverse_order = n : -1 : 1;
Fac.Q = Fac.Q(:, reverse_order);
R = R(reverse_order, :);
Fac.re_size = n - s;
Fac.col_sk_size = s;
Fac.row_sk_size = m - Fac.re_size;
Fac.V_sk = R((Fac.re_size + 1) : n, :);
data.Amat = data.Amat * Fac.Q;

% Eliminate redundant columns and retain the triangular solve blocks.
[Fac.P, R] = qr(data.Amat(:, 1 : Fac.re_size));
Fac.A_re_re = R(1 : Fac.re_size, :);
A_sk = Fac.P' * data.Amat(:, (Fac.re_size + 1) : n);
Fac.A_re_sk = A_sk(1 : Fac.re_size, :);
data.Amat = A_sk((Fac.re_size + 1) : m, :);

% Retain U_re for the solve and pass the remaining generators upward.
U = Fac.P' * data.Umat;
Fac.U_re = U(1 : Fac.re_size, :);
data.Umat = U((Fac.re_size + 1) : m, :);
data.Vmat = Fac.V_sk;

end
