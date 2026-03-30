function ConstructRootGenerators2(A, xy, col_pos)

arguments (Input)
    A NUDFT2_2D_HSS;
    xy (:, 2) double;
    col_pos (:, 2) double;
end

M = A.row_global_size_;
nx = A.nx_;
ny = A.ny_;
N = A.col_global_size_;
x_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, nx);
y_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, ny);

if A.leaf_ == 1
    A.row_ind_ = A.row_offset_ + (1 : A.row_size_)';
    A.col_ind_ = A.col_offset_ + (1 : A.col_size_)';
    
    xy_I = xy(A.row_ind_, :);
    gamma_x_I = exp(-2 * pi * 1i * xy_I(:, 1));
    gamma_y_I = exp(-2 * pi * 1i * xy_I(:, 2));
    col_pos_J = col_pos(A.col_ind_, :);
    xi_x_J = exp(-2 * pi * 1i * col_pos_J(:, 1) / nx);
    xi_y_J = exp(-2 * pi * 1i * col_pos_J(:, 2) / ny);

    % Assign full mat.
    Ax = x_kernel_fun(gamma_x_I, xi_x_J);
    Ay = y_kernel_fun(gamma_y_I, xi_y_J);
    A.Amat_ = Ax .* Ay;

    % Clear.
    A.row_ind_ = [];
    A.col_ind_ = [];
else
    % Construct B.
    for i = 1 : A.num_children_
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

    % Clear data.
    for i = 1 : A.num_children_
        A.children_{i}.row_ind_ = [];
        A.children_{i}.col_ind_ = [];
    end
end

end