function URV_RootFactor(A)
% URV_RootFactor

% Jingyu Liu, November 22, 2024.

arguments (Input)
    A HSS;
end

if A.leaf_ == 0
    A.URV_Merge();
end

m = size(A.Amat_, 1);
n = size(A.Amat_, 2);
% fprintf("m: %d, n: %d\n", m, n);

[A.URV_P_, A.URV_A_re_re_] = qr(A.Amat_);

% Clear data.
A.Amat_ = [];

end