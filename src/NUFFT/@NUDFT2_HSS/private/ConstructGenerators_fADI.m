function data = ConstructGenerators_fADI(A, x, tol)
% Construct generators recursively using points in global leaf order.

N = A.col_global_size_;
pos_start = A.col_offset_;
pos_end = pos_start + A.col_size_ - 1;
data = struct;

if A.leaf_ == 1
    data.row_x = x(A.row_offset_ + (1 : A.row_size_));
    data.col_pos = (pos_start : pos_end)';
else
    children = cell(1, A.num_children_);
    row_size = 0;
    col_size = 0;
    for i = 1 : A.num_children_
        children{i} = ConstructGenerators_fADI(A.children_{i}, x, tol);
        row_size = row_size + A.children_{i}.row_rank_;
        col_size = col_size + A.children_{i}.col_rank_;
    end

    A.Bmat_ = cell(A.num_children_);
    for i = 1 : A.num_children_
        z_I = exp(-2 * pi * 1i * children{i}.row_x);
        for j = 1 : A.num_children_
            if i ~= j
                w_J = exp(-2 * pi * 1i * children{j}.col_pos / N);
                A.Bmat_{i, j} = NUFFT2_Kernel(z_I, w_J, N);
            end
        end
    end
    if A.level_ == 0
        return;
    end

    data.row_x = zeros(row_size, 1);
    data.col_pos = zeros(col_size, 1);
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        data.row_x(row_offset + (1 : A.children_{i}.row_rank_)) = children{i}.row_x;
        data.col_pos(col_offset + (1 : A.children_{i}.col_rank_)) = children{i}.col_pos;
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end
end

z_I = exp(-2 * pi * 1i * data.row_x);
w_J = exp(-2 * pi * 1i * data.col_pos / N);
if A.level_ == 0
    A.Amat_ = NUFFT2_Kernel(z_I, w_J, N);
    return;
end

% Construct the row interpolation basis from fADI factors.
x = data.row_x;
tmp_ind = find(x >= (1 - 1 / N / 2));
x(tmp_ind) = x(tmp_ind) - 1;
Ir = [...
    exp(-2 * pi * 1i * min(x)), ...
    exp(-2 * pi * 1i * max(x)), ...
    exp(-2 * pi * 1i * (pos_end + 1) / N), ...
    exp(-2 * pi * 1i * (pos_start - 1) / N)];
[~, ~, ~, cr] = mobiusT(Ir);
k = ceil(1/pi^2*log(4/tol)*log(16*cr));
[alpha, beta] = getshifts_adi(Ir, k);
[row_sk, U, A.row_rank_] = fADI_Row_NUDFT2(z_I, z_I.^N - 1, k, alpha, beta);

% Construct the column interpolation basis from fADI factors.
Ic = [...
    exp(-2 * pi * 1i * (pos_end + 1 / 2) / N), ...
    exp(-2 * pi * 1i * (pos_start - 1 / 2) / N), ...
    exp(-2 * pi * 1i * min(data.col_pos) / N), ...
    exp(-2 * pi * 1i * max(data.col_pos) / N)];
[~, ~, ~, cc] = mobiusT(Ic);
s = ceil(1/pi^2*log(4/tol)*log(16*cc));
[alpha, beta] = getshifts_adi(Ic, s);
[col_sk, V, A.col_rank_] = fADI_Col_NUDFT2(w_J, conj(w_J) / N, s, alpha, beta);

if A.leaf_ == 1
    A.Umat_ = U;
    A.Vmat_ = V;
    A.Amat_ = NUFFT2_Kernel(z_I, w_J, N);
else
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        A.children_{i}.Rmat_ = U(row_offset + (1 : A.children_{i}.row_rank_), :);
        A.children_{i}.Wmat_ = V(col_offset + (1 : A.children_{i}.col_rank_), :);
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end
end
data.row_x = data.row_x(row_sk);
data.col_pos = data.col_pos(col_sk);

end
