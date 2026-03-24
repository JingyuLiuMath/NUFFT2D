function ConstructRootGenerators(A)
% ConstructGenerators

% Jingyu Liu, November 20, 2024.

arguments (Input)
    A NUDFT2_HSS;
end

N = A.col_global_size_;
kernel_fun = @(z, w) NUFFT2_Kernel(z, w, N);

if A.leaf_ == 1
    z_I = exp(-2 * pi * 1i * A.row_x_);
    w_J = exp(-2 * pi * 1i * A.col_pos_ / N);

    % Assign full mat.
    A.Amat_ = kernel_fun(z_I, w_J);

    % Clear.
    A.row_x_ = [];
    A.col_pos_ = [];
else
    % Construct B.
    for i = 1 : A.num_children_
        c_z_I = exp(-2 * pi * 1i * A.children_{i}.row_x_);
        for j = [1 : (i - 1), (i + 1) : A.num_children_]
            c_w_J = exp(-2 * pi * 1i * A.children_{j}.col_pos_ / N);
            A.Bmat_{i, j} = kernel_fun(c_z_I, c_w_J);
        end
    end

    % Clear data.
    for i = 1 : A.num_children_
        A.children_{i}.row_x_ = [];
        A.children_{i}.col_pos_ = [];
    end
end

end