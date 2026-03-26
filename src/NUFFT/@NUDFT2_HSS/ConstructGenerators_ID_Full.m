function ConstructGenerators_ID_Full(A, ...
    x, col_pos, ...
    row_inact, col_inact, ...
    level, rank_or_tol)

arguments (Input)
    A NUDFT2_HSS;
    x (:, 1) double;
    col_pos (:, 1) double;
    row_inact (:, 1) double;
    col_inact (:, 1) double;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

M = A.row_global_size_;
N = A.col_global_size_;
kernel_fun = @(z, w) NUFFT2_Kernel(z, w, N);

if A.level_ == level
    if A.leaf_ == 1
        A.row_ind_ = A.row_offset_ + (1 : A.row_size_)';
        A.col_ind_ = A.col_offset_ + (1 : A.col_size_)';
    else
        % Merge data from children and construct B.
        A.row_ind_ = [];
        A.col_ind_ = [];
        for i = 1 : A.num_children_
            A.row_ind_ = [A.row_ind_; A.children_{i}.row_ind_];
            A.col_ind_ = [A.col_ind_; A.children_{i}.col_ind_];
            c_x_I = x(A.children_{i}.row_ind_, :);
            c_z_I = exp(-2 * pi * 1i * c_x_I);
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                c_col_pos_J = col_pos(A.children_{j}.col_ind_, :);
                c_w_J = exp(-2 * pi * 1i * c_col_pos_J / N);
                A.Bmat_{i, j} = kernel_fun(c_z_I, c_w_J);
            end
        end
    end

    % Construct U.
    x_I = x(A.row_ind_, :);
    z_I = exp(-2 * pi * 1i * x_I);
    Jc = (1 : N)';
    self_col_ind = A.col_offset_ + (1 : A.col_size_)';
    Jc([self_col_ind; col_inact]) = [];
    col_pos_Jc = col_pos(Jc, :);
    w_Jc = exp(-2 * pi * 1i * col_pos_Jc / N);
    A_I_Jc = kernel_fun(z_I, w_Jc);
    [row_sk, U, A.row_rank_, row_re] = LowRank_Row_ID(A_I_Jc, rank_or_tol);

    % Construct V using proxy surface.
    col_pos_J = col_pos(A.col_ind_, :);
    w_J = exp(-2 * pi * 1i * col_pos_J / N);
    Ic = (1 : M)';
    self_row_ind = A.row_offset_ + (1 : A.row_size_)';
    Ic([self_row_ind; row_inact]) = [];
    x_Ic = x(Ic, :);
    z_Ic = exp(-2 * pi * 1i * x_Ic);
    A_Ic_J = kernel_fun(z_Ic, w_J);
    [col_sk, V, A.col_rank_, col_re] = LowRank_ID(A_Ic_J, rank_or_tol);

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
            current_size = size(A.children_{i}.row_ind_, 1);
            A.Rmat_{i} = U((offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end

        offset = 0;
        for i = 1 : A.num_children_
            current_size = size(A.children_{i}.col_ind_, 1);
            A.Wmat_{i} = V((offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end

        % Clear data.
        for i = 1 : A.num_children_
            A.children_{i}.row_ind_ = [];
            A.children_{i}.col_ind_ = [];
        end
    end

    % Update row and col.
    A.row_re_ = A.row_ind_(row_re);
    A.col_re_ = A.col_ind_(col_re);
    A.row_ind_ = A.row_ind_(row_sk);
    A.col_ind_ = A.col_ind_(col_sk);
else
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_ID_Full(...
            x, col_pos, ...
            row_inact, col_inact, ...
            level, rank_or_tol);
    end
end

end