function ConstructGenerators_fADI(A, level, rank_or_tol)

arguments (Input)
    A NUDFT2_HSS;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

N = A.col_global_size_;
kernel_fun = @(z, w) NUFFT2_Kernel(z, w, N);

if A.level_ == level
    if A.leaf_ == 0
        % Merge data from children and construct B.
        A.row_x_ = [];
        A.col_pos_ = [];
        for i = 1 : A.num_children_
            A.row_x_ = [A.row_x_; A.children_{i}.row_x_];
            A.col_pos_ = [A.col_pos_; A.children_{i}.col_pos_];
            c_z_I = exp(-2 * pi * 1i * A.children_{i}.row_x_);
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                c_w_J = exp(-2 * pi * 1i * A.children_{j}.col_pos_ / N);
                A.Bmat_{i, j} = kernel_fun(c_z_I, c_w_J);
            end
        end
    end

    % Construct U or R.
    z_I = exp(-2 * pi * 1i * A.row_x_);  % row points.
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
    if rank_or_tol >= 1
        k = min(rank_or_tol, A.row_size_);
    else
        k = ceil(1/pi^2*log(4/rank_or_tol)*log(16*cr));
    end
    [alpha, beta] = getshifts_adi(Ir, k);
    [row_sk, U, A.row_rank_] = fADI_Row_NUDFT2(...
        z_I, u, k, alpha, beta);

    % A_I_Jc = kernel_fun_real(A.row_x_, A.real_res_col_pos_);
    % [row_sk, U, A.row_rank_] = LowRank_Row_ID(A_I_Jc, k);

    % Construct V or W.
    w_J = exp(-2 * pi * 1i * A.col_pos_ / N);
    v = conj(w_J) / N;
    Ic = [...
        exp(-2 * pi * 1i * (A.pos_end_ + 1 / 2) / N), ...
        exp(-2 * pi * 1i * (A.pos_start_ - 1 / 2) / N), ...
        exp(-2 * pi * 1i * min(A.col_pos_) / N), ...
        exp(-2 * pi * 1i * max(A.col_pos_) / N)];
    [~, ~, ~, cc] = mobiusT(Ic);
    if rank_or_tol >= 1
        s = min(rank_or_tol, A.col_size_);
    else
        s = ceil(1/pi^2*log(4/rank_or_tol)*log(16*cc));
    end
    [alpha, beta] = getshifts_adi(Ic, s);
    [col_sk, V, A.col_rank_] = fADI_Col_NUDFT2(...
        w_J, v, s, alpha, beta);

    % A_Ic_J = kernel_fun_real(A.real_res_row_x_, A.col_pos_);
    % [col_sk, V, A.col_rank_] = LowRank_ID(A_Ic_J, s);

    if A.leaf_ == 1
        % Assign U and V.
        A.Umat_ = U;
        A.Vmat_ = V;

        % Assign full mat.
        A.Amat_ = kernel_fun(z_I, w_J);
    else
        % Assign R and W.
        offset = 0;
        for i = 1 : A.num_children_
            current_size = size(A.children_{i}.row_x_, 1);
            A.Rmat_{i} = U((offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end
        
        offset = 0;
        for i = 1 : A.num_children_
            current_size = size(A.children_{i}.col_pos_, 1);
            A.Wmat_{i} = V((offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
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

else
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_fADI(level, rank_or_tol);
    end
end

end