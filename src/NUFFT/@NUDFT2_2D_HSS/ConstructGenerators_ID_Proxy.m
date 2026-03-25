function ConstructGenerators_ID_Proxy(A, level, rank_or_tol)
% ConstructGenerators_ID_Proxy

% Jingyu Liu, February 4, 2025.

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

M = A.row_global_size_;
nx = A.nx_;
ny = A.ny_;
N = A.col_global_size_;
x_half_length = 0.5 / nx;
y_half_length = 0.5 / ny;
x_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, nx);
y_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, ny);

if A.level_ == level
    if A.leaf_ == 0
        A.row_xy_ = [];
        A.col_pos_ = [];
        for i = 1 : A.num_children_
            A.row_xy_ = [A.row_xy_; A.children_{i}.row_xy_];
            A.col_pos_ = [A.col_pos_; A.children_{i}.col_pos_];
            c_gamma_x_I = exp(-2 * pi * 1i * A.children_{i}.row_xy_(:, 1));
            c_gamma_y_I = exp(-2 * pi * 1i * A.children_{i}.row_xy_(:, 2));
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                c_xi_x_J = exp(-2 * pi * 1i * A.children_{j}.col_pos_(:, 1) / nx);
                c_xi_y_J = exp(-2 * pi * 1i * A.children_{j}.col_pos_(:, 2) / ny);
                Ax = x_kernel_fun(c_gamma_x_I, c_xi_x_J);
                Ay = y_kernel_fun(c_gamma_y_I, c_xi_y_J);
                A.Bmat_{i, j} = Ax .* Ay;
            end
        end
    end

    nx_ny = max(nx, ny);
    proxy_size = 5;

    % Generate proxy surface.
    % x_all = (A.x_pos_start_ : A.x_pos_end_)' / nx;
    % x_all = x_all + 1e-12 * randn(size(x_all));
    x_all = linspace(...
        A.x_pos_start_ / nx - x_half_length, ...
        A.x_pos_end_ / nx + x_half_length, ...
        A.x_col_size_);
    tmp_ind = find(x_all >= 1);
    x_all(tmp_ind) = x_all(tmp_ind) - 1;
    tmp_ind = find(x_all < 0);
    x_all(tmp_ind) = x_all(tmp_ind) + 1;

    col_pt_x_start = A.x_pos_start_ / nx;
    col_pt_x_end = A.x_pos_end_ / nx;
    x_proxy = GenerateProxySurface(nx, ...
        col_pt_x_start - x_half_length, ...
        col_pt_x_end + x_half_length, ...
        proxy_size);
    tmp_ind = find(x_proxy >= 1);
    x_proxy(tmp_ind) = x_proxy(tmp_ind) - 1;
    tmp_ind = find(x_proxy < 0);
    x_proxy(tmp_ind) = x_proxy(tmp_ind) + 1;

    % y_all = (A.y_pos_start_ : A.y_pos_end_)' / ny;
    % y_all = y_all + 1e-12 * randn(size(y_all));
    y_all = linspace(...
        A.y_pos_start_ / ny - y_half_length, ...
        A.y_pos_end_ / ny + y_half_length, ...
        A.y_col_size_);
    tmp_ind = find(y_all >= 1);
    y_all(tmp_ind) = y_all(tmp_ind) - 1;
    tmp_ind = find(y_all < 0);
    y_all(tmp_ind) = y_all(tmp_ind) + 1;

    col_pt_y_start = A.y_pos_start_ / ny;
    col_pt_y_end = A.y_pos_end_ / ny;
    y_proxy = GenerateProxySurface(ny, ...
        col_pt_y_start - y_half_length, ...
        col_pt_y_end + y_half_length, ...
        proxy_size);
    tmp_ind = find(y_proxy >= 1);
    y_proxy(tmp_ind) = y_proxy(tmp_ind) - 1;
    tmp_ind = find(y_proxy < 0);
    y_proxy(tmp_ind) = y_proxy(tmp_ind) + 1;

    proxy_surface_lr = TensorProduct2D(x_proxy, y_all);
    proxy_surface_tb = TensorProduct2D(x_all, y_proxy);
    proxy_surface_real = [proxy_surface_lr; proxy_surface_tb];
    proxy_surface = [...
        exp(-2 * pi * 1i * proxy_surface_real(:, 1)), ...
        exp(-2 * pi * 1i * proxy_surface_real(:, 2))];

    figure();
    plot(A.row_xy_(:, 1), A.row_xy_(:, 2), "rx", "DisplayName", "row pts");
    hold on;
    plot(A.col_pos_(:, 1) / nx, A.col_pos_(:, 2) / ny, "go", "DisplayName", "col pts");
    plot(proxy_surface_real(:, 1), proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
    legend;
    keyboard;

    % Construct U using proxy surface.
    gamma_x_I = exp(-2 * pi * 1i * A.row_xy_(:, 1));
    gamma_y_I = exp(-2 * pi * 1i * A.row_xy_(:, 2));
    Ax_I_proxy = x_kernel_fun(gamma_x_I, proxy_surface(:, 1));
    Ay_I_proxy = y_kernel_fun(gamma_y_I, proxy_surface(:, 2));
    A_I_proxy = Ax_I_proxy .* Ay_I_proxy;
    [row_sk, U, A.row_rank_] ...
        = LowRank_Row_ID(A_I_proxy, rank_or_tol);

    % Construct V using proxy surface.
    xi_x_J = exp(-2 * pi * 1i * A.col_pos_(:, 1) / nx);
    xi_y_J = exp(-2 * pi * 1i * A.col_pos_(:, 2) / ny);
    Ax_proxy_J = x_kernel_fun(proxy_surface(:, 1), xi_x_J);
    Ay_proxy_J = y_kernel_fun(proxy_surface(:, 2), xi_y_J);
    A_proxy_J = Ax_proxy_J .* Ay_proxy_J;
    [col_sk, V, A.col_rank_] ...
        = LowRank_ID(A_proxy_J, rank_or_tol);

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
            current_size = size(A.children_{i}.row_xy_, 1);
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
            A.children_{i}.row_xy_ = [];
            A.children_{i}.col_pos_ = [];
        end
    end

    % Update row and col.
    A.row_xy_ = A.row_xy_(row_sk, :);
    A.col_pos_ = A.col_pos_(col_sk, :);
else
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_ID_Proxy(level, rank_or_tol);
    end
end

end