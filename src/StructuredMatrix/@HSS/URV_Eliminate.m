function URV_Eliminate(A, level)
% URV_Eliminate

% Jingyu Liu, November 22, 2024.

arguments (Input)
    A HSS;
    level (1, 1) double;
end

A.urv_leaf_ = A.leaf_;

if A.level_ == level
    if A.leaf_ == 0
        A.URV_Merge();
    end
    
    m = size(A.Umat_, 1);
    n = size(A.Vmat_, 1);
    k = A.row_rank_;
    s = A.col_rank_;

    % Size reduction.
    if m >= (k + n)
        [A.URV_Omega_, R] = qr([A.Umat_, A.Amat_], "econ");
        A.URV_m_old_ = m;
        m = k + n;
        A.URV_m_new_ = m;
        A.Umat_ = R(1 : m, 1 : k);
        A.Amat_ = R(1 : m, (k + 1) : (k + n));
    end

    % Elimination.
    % *********************************************************************
    % Step 1. Zeroing out V.
    [A.URV_Q_, R] = qr(A.Vmat_);  % V: n-by-s, Q: n-by-n, R: n-by-s.
    reverse_order = n : -1 : 1;
    A.URV_Q_ = A.URV_Q_(:, reverse_order);  % Q: n-by-n.
    R = R(reverse_order, :);
    re_size = n - s;  % m >= n - s.
    col_sk_size = n - re_size;
    row_sk_size = m - re_size;
    A.URV_re_size_ = re_size;
    A.URV_col_sk_size_ = col_sk_size;
    A.URV_row_sk_size_ = row_sk_size;
    A.URV_V_sk_ = R((re_size + 1) : n, :);  % V_sk: s-by-s.
    % Update A.
    A.Amat_ = A.Amat_ * A.URV_Q_;
    % ---------------------------------------------------------------------

    % *********************************************************************
    % Step 2. QR on Amat.
    A_all_1 = A.Amat_(:, 1 : re_size);  % A_all_1: m-by-(n - s).
    A_all_2 = A.Amat_(:, (re_size + 1) : n);  % A_all_2: m-by-s.
    [A.URV_P_, R] = qr(A_all_1);  % P: m-by-m, R: m -by-(n - s).
    A.URV_A_re_re_ = R(1 : re_size, :);  % A_re_re_: (n - s)-by-(n - s).

    A_all_2 = A.URV_P_' * A_all_2;
    A.URV_A_re_sk_ ...
        = A_all_2(1 : re_size, :);  % A_re_sk: (n - s)-by-s.
    A.URV_A_sk_sk_ ...
        = A_all_2((re_size + 1) : m, :);  % A_sk_sk: (m - n + s)-by-s.
    % ---------------------------------------------------------------------

    % *********************************************************************
    % Step 3. Update U.
    A.Umat_ = A.URV_P_' * A.Umat_;
    A.URV_U_re_ = A.Umat_(1 : re_size, :);  % U_re: (n - s)-by-k.
    A.URV_U_sk_ = A.Umat_((re_size + 1) : m, :);  % U_sk: (m - n + s)-by-k.
    % ---------------------------------------------------------------------

    % *********************************************************************
    % Clear data.
    A.Amat_ = [];
    A.Umat_ = [];
    A.Vmat_ = [];
    % ---------------------------------------------------------------------
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.URV_Eliminate(level);
    end
end

end