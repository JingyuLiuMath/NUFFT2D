function ConstructRootGenerators2(A, x, col_pos)

arguments (Input)
    A NUDFT2_HSS;
    x (:, 1) double;
    col_pos (:, 1) double;
end

N = A.col_global_size_;
kernel_fun = @(z, w) NUFFT2_Kernel(z, w, N);

if A.leaf_ == 1
    A.row_ind_ = A.row_offset_ + (1 : A.row_size_)';
    A.col_ind_ = A.col_offset_ + (1 : A.col_size_)';

    x_I = x(A.row_ind_, :);
    z_I = exp(-2 * pi * 1i * x_I);
    col_pos_J = col_pos(A.col_ind_, :);
    w_J = exp(-2 * pi * 1i * col_pos_J / N);

    % Assign full mat.
    A.Amat_ = kernel_fun(z_I, w_J);

    % Clear.
    A.row_ind_ = [];
    A.col_ind_ = [];
else
    % Construct B.
    for i = 1 : A.num_children_
        c_x_I = x(A.children_{i}.row_ind_, :);
        c_z_I = exp(-2 * pi * 1i * c_x_I);
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            c_col_pos_J = col_pos(A.children_{j}.col_ind_, :);
            c_w_J = exp(-2 * pi * 1i * c_col_pos_J / N);
            A.Bmat_{i, j} = kernel_fun(c_z_I, c_w_J);
        end
    end

    % Clear data.
    for i = 1 : A.num_children_
        A.children_{i}.row_ind_ = [];
        A.children_{i}.col_ind_ = [];
    end
end

end