function ConstructGenerators_ID_Full(A, ...
    xy, col_pos, ...
    row_inact, col_inact, ...
    level, rank_or_tol)

arguments (Input)
    A NUDFT2_2D_HSS;
    xy (:, 2) double;
    col_pos (:, 2) double;
    row_inact (:, 1) double;
    col_inact (:, 1) double;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

M = A.row_global_size_;
nx = A.nx_;
ny = A.ny_;
N = A.col_global_size_;
x_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, nx);
y_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, ny);

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
            c_xy_I = xy(A.children_{i}.row_ind_, :);
            c_gamma_x_I = exp(-2 * pi * 1i * c_xy_I(:, 1));
            c_gamma_y_I = exp(-2 * pi * 1i * c_xy_I(:, 2));
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                c_col_pos_J = col_pos(A.children_{j}.col_ind_, :);
                c_xi_x_J = exp(-2 * pi * 1i * c_col_pos_J(:, 1) / nx);
                c_xi_y_J = exp(-2 * pi * 1i * c_col_pos_J(:, 2) / ny);
                Ax = x_kernel_fun(c_gamma_x_I, c_xi_x_J);
                Ay = y_kernel_fun(c_gamma_y_I, c_xi_y_J);
                A.Bmat_{i, j} = Ax .* Ay;
            end
        end
    end

    % Construct U.
    xy_I = xy(A.row_ind_, :);
    gamma_x_I = exp(-2 * pi * 1i * xy_I(:, 1));
    gamma_y_I = exp(-2 * pi * 1i * xy_I(:, 2));
    Jc = (1 : N)';
    self_col_ind = A.col_offset_ + (1 : A.col_size_)';
    Jc([self_col_ind; col_inact]) = [];
    col_pos_Jc = col_pos(Jc, :);
    xi_x_Jc = exp(-2 * pi * 1i * col_pos_Jc(:, 1) / nx);
    xi_y_Jc = exp(-2 * pi * 1i * col_pos_Jc(:, 2) / ny);
    Ax_I_Jc = x_kernel_fun(gamma_x_I, xi_x_Jc);
    Ay_I_Jc = y_kernel_fun(gamma_y_I, xi_y_Jc);
    A_I_Jc = Ax_I_Jc .* Ay_I_Jc;
    [row_sk, U, A.row_rank_, row_re] = LowRank_Row_ID(A_I_Jc, rank_or_tol);

    % Construct V.
    col_pos_J = col_pos(A.col_ind_, :);
    xi_x_J = exp(-2 * pi * 1i * col_pos_J(:, 1) / nx);
    xi_y_J = exp(-2 * pi * 1i * col_pos_J(:, 2) / ny);
    Ic = (1 : M)';
    self_row_ind = A.row_offset_ + (1 : A.row_size_)';
    Ic([self_row_ind; row_inact]) = [];
    xy_Ic = xy(Ic, :);
    gamma_x_Ic = exp(-2 * pi * 1i * xy_Ic(:, 1));
    gamma_y_Ic = exp(-2 * pi * 1i * xy_Ic(:, 2));
    Ax_Ic_J = x_kernel_fun(gamma_x_Ic, xi_x_J);
    Ay_Ic_J = y_kernel_fun(gamma_y_Ic, xi_y_J);
    A_Ic_J = Ax_Ic_J .* Ay_Ic_J;
    [col_sk, V, A.col_rank_, col_re] = LowRank_ID(A_Ic_J, rank_or_tol);

    % figure();
    % plot(xy_I(:, 1), xy_I(:, 2), "rx", "DisplayName", "row pts");
    % hold on;
    % plot(col_pos_J(:, 1) / nx, col_pos_J(:, 2) / ny, "go", "DisplayName", "col pts");
    % keyboard;

    if A.leaf_ == 1
        % Assign U and V.
        A.Umat_ = U;
        A.Vmat_ = V;

        % Assign full mat.
        Ax = x_kernel_fun(gamma_x_I, xi_x_J);
        Ay = y_kernel_fun(gamma_y_I, xi_y_J);
        A.Amat_ = Ax .* Ay;
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
            xy, col_pos, ...
            row_inact, col_inact, ...
            level, rank_or_tol);
    end
end

end