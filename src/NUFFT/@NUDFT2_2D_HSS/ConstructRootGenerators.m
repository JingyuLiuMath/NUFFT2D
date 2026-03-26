function ConstructRootGenerators(A)

arguments (Input)
    A NUDFT2_2D_HSS;
end

nx = A.nx_;
ny = A.ny_;
x_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, nx);
y_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, ny);
x_kernel_fun_real = @(x, pos) NUFFT2_Kernel_Real(x, pos, nx);
y_kernel_fun_real = @(x, pos) NUFFT2_Kernel_Real(x, pos, ny);

if A.leaf_ == 1
    % z_Ix = exp(-2 * pi * 1i * A.row_xy_(:, 1));
    % z_Iy = exp(-2 * pi * 1i * A.row_xy_(:, 2));
    % w_Jx = exp(-2 * pi * 1i * A.col_pos_(:, 1) / nx);
    % w_Jy = exp(-2 * pi * 1i * A.col_pos_(:, 2) / ny);

    % Assign full mat.
    % Ax = x_kernel_fun(z_Ix, w_Jx);
    % Ay = y_kernel_fun(z_Iy, w_Jy);
    % A.Amat_ = Ax .* Ay;
    Ax = x_kernel_fun_real(A.row_xy_(:, 1), A.col_pos_(:, 1) / nx);
    Ay = y_kernel_fun_real(A.row_xy_(:, 2), A.col_pos_(:, 2) / ny);
    A.Amat_ = Ax .* Ay;
else
    % Construct B.
    for i = 1 : A.num_children_
        % c_z_Ix = exp(-2 * pi * 1i * A.children_{i}.row_x_(:, 1));
        % c_z_Iy = exp(-2 * pi * 1i * A.children_{i}.row_x_(:, 2));
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            % c_w_Jx = exp(-2 * pi * 1i * A.children_{j}.col_pos_(:, 1) / nx);
            % c_w_Jy = exp(-2 * pi * 1i * A.children_{j}.col_pos_(:, 2) / ny);
            % Ax = x_kernel_fun(c_z_Ix, c_w_Jx);
            % Ay = y_kernel_fun(c_z_Iy, c_w_Jy);
            % A.Bmat_{i, j} = Ax .* Ay;
            Ax = x_kernel_fun_real(...
                A.children_{i}.row_xy_(:, 1), ...
                A.children_{j}.col_pos_(:, 1) / nx);
            Ay = y_kernel_fun_real(...
                A.children_{i}.row_xy_(:, 2), ...
                A.children_{j}.col_pos_(:, 2) / ny);
            A.Bmat_{i, j} = Ax .* Ay;
        end
    end

    % Clear data.
    for i = 1 : A.num_children_
        A.children_{i}.row_xy_ = [];
        A.children_{i}.col_pos_ = [];
    end
end

end