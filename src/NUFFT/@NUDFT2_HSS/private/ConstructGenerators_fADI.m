function ConstructGenerators_fADI(A, tol)

arguments (Input)
    A NUDFT2_HSS;
    tol (1, 1) double;
end

M = A.row_global_size_;
N = A.col_global_size_;
kernel_fun = @(z, w) NUFFT2_Kernel(z, w, N);

% Construct children first.
for i = 1 : A.num_children_
    ConstructGenerators_fADI(A.children_{i}, tol);
end

if A.level_ == 0
    ConstructRootGenerators(A);
    return;
end

if A.leaf_ == 0
    % Merge data from children and construct B.
    row_size = 0;
    col_size = 0;
    for i = 1 : A.num_children_
        row_size = row_size + A.children_{i}.row_rank_;
        col_size = col_size + A.children_{i}.col_rank_;
    end
    A.row_x_ = zeros(row_size, 1);
    A.col_pos_ = zeros(col_size, 1);
    A.Bmat_ = cell(A.num_children_);
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        A.row_x_(row_offset + (1 : A.children_{i}.row_rank_)) ...
            = A.children_{i}.row_x_;
        A.col_pos_(col_offset + (1 : A.children_{i}.col_rank_)) ...
            = A.children_{i}.col_pos_;
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;

        c_z_I = exp(-2 * pi * 1i * A.children_{i}.row_x_);
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            c_w_J = exp(-2 * pi * 1i * A.children_{j}.col_pos_ / N);
            A.Bmat_{i, j} = NUFFT2_Kernel(c_z_I, c_w_J, N);
        end
    end
end

% Construct U or R.
z_I = exp(-2 * pi * 1i * A.row_x_);
u = z_I.^N - 1;
row_x = A.row_x_;
tmp_ind = find(row_x >= (1 - 1 / N / 2));
row_x(tmp_ind) = row_x(tmp_ind) - 1;
Ir = [...
    exp(-2 * pi * 1i * min(row_x)), ...
    exp(-2 * pi * 1i * max(row_x)), ...
    exp(-2 * pi * 1i * (A.pos_end_ + 1) / N), ...
    exp(-2 * pi * 1i * (A.pos_start_ - 1) / N)];
[~, ~, ~, cr] = mobiusT(Ir);
k = ceil(1/pi^2*log(4/tol)*log(16*cr));
[alpha, beta] = getshifts_adi(Ir, k);
[row_sk, U, A.row_rank_] = fADI_Row_NUDFT2(z_I, u, k, alpha, beta);

% Construct V or W.
w_J = exp(-2 * pi * 1i * A.col_pos_ / N);
v = conj(w_J) / N;
Ic = [...
    exp(-2 * pi * 1i * (A.pos_end_ + 1 / 2) / N), ...
    exp(-2 * pi * 1i * (A.pos_start_ - 1 / 2) / N), ...
    exp(-2 * pi * 1i * min(A.col_pos_) / N), ...
    exp(-2 * pi * 1i * max(A.col_pos_) / N)];
[~, ~, ~, cc] = mobiusT(Ic);
s = ceil(1/pi^2*log(4/tol)*log(16*cc));
[alpha, beta] = getshifts_adi(Ic, s);
[col_sk, V, A.col_rank_] = fADI_Col_NUDFT2(w_J, v, s, alpha, beta);

if A.leaf_ == 1
    % Assign U and V.
    A.Umat_ = U;
    A.Vmat_ = V;

    % Assign full mat.
    A.Amat_ = kernel_fun(z_I, w_J);
else
    % Assign R and W.
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        A.children_{i}.Rmat_ = U(row_offset + (1 : A.children_{i}.row_rank_), :);
        A.children_{i}.Wmat_ = V(col_offset + (1 : A.children_{i}.col_rank_), :);
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end

    % Clear data.
    for i = 1 : A.num_children_
        A.children_{i}.row_x_ = [];
        A.children_{i}.col_pos_ = [];
    end
end

% Update row and col.
A.row_x_ = A.row_x_(row_sk);
A.col_pos_ = A.col_pos_(col_sk);

end
