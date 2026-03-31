function ConstructGenerators_ID_Proxy(A, level, rank_or_tol)

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

M = A.row_global_size_;
nx = A.nx_;
ny = A.ny_;
N = A.col_global_size_;
x_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, nx);
y_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, ny);
hx = 1 / nx;
hy = 1 / ny;
nx_ny = max(nx, ny);
proxy_layer_size = log2(nx_ny);
sampling_size = 2 * nx_ny * proxy_layer_size;

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

    % Generate proxy surface.
    row_proxy_surface_real = [];
    col_proxy_surface_real = [];

    x_pt_start = A.x_pos_start_ / nx;
    x_pt_end = A.x_pos_end_ / nx ;
    y_pt_start = A.y_pos_start_ / ny;
    y_pt_end = A.y_pos_end_ / ny;

    x_self_proxy_size = A.x_col_size_;
    x_self_pts = linspace(x_pt_start, x_pt_end, x_self_proxy_size)';
    x_self_pts = [...
        x_self_pts -  proxy_layer_size * hx;
        x_self_pts;
        x_self_pts +  proxy_layer_size * hx];

    y_self_proxy_size = A.y_col_size_;
    y_self_pts = linspace(y_pt_start, y_pt_end, y_self_proxy_size)';
    y_self_pts = [...
        y_self_pts -  proxy_layer_size * hy;
        y_self_pts;
        y_self_pts +  proxy_layer_size * hy];
    
    for ity = 1 : proxy_layer_size
        row_proxy_surface_real = [row_proxy_surface_real; ...
            [x_self_pts, (y_pt_start - ity * hy) * ones(size(x_self_pts))]];
        row_proxy_surface_real = [row_proxy_surface_real; ...
            [x_self_pts, (y_pt_end + ity * hy) * ones(size(x_self_pts))]];
    end
    for itx = 1 : proxy_layer_size
        row_proxy_surface_real = [row_proxy_surface_real; ...
            [(x_pt_start - itx * hx) * ones(size(y_self_pts)), y_self_pts]];
        row_proxy_surface_real = [row_proxy_surface_real; ...
            [(x_pt_end + itx * hx) * ones(size(y_self_pts)), y_self_pts]];
    end

    tmp_ind = row_proxy_surface_real(:, 1) < 0;
    row_proxy_surface_real(tmp_ind, 1) = row_proxy_surface_real(tmp_ind, 1) + 1;
    tmp_ind = row_proxy_surface_real(:, 1) >= 1;
    row_proxy_surface_real(tmp_ind, 1) = row_proxy_surface_real(tmp_ind, 1) - 1;

    tmp_ind = row_proxy_surface_real(:, 2) < 0;
    row_proxy_surface_real(tmp_ind, 2) = row_proxy_surface_real(tmp_ind, 2) + 1;
    tmp_ind = row_proxy_surface_real(:, 2) >= 1;
    row_proxy_surface_real(tmp_ind, 2) = row_proxy_surface_real(tmp_ind, 2) - 1;
    
    row_proxy_surface = [...
        exp(-2 * pi * 1i * row_proxy_surface_real(:, 1)), ...
        exp(-2 * pi * 1i * row_proxy_surface_real(:, 2))];

    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(sampling_size, ...
        x_pt_start, x_pt_end, ...
        y_pt_start - proxy_layer_size * hy, y_pt_start - hy)];
    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(sampling_size, ...
        x_pt_start, x_pt_end, ...
        y_pt_end + hy, y_pt_end + proxy_layer_size * hy)];
    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(sampling_size, ...
        x_pt_start - proxy_layer_size * hx, x_pt_start - hx, ...
        y_pt_start, y_pt_end)];
    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(sampling_size, ...
        x_pt_end + hx, x_pt_end + proxy_layer_size * hx, ...
        y_pt_start, y_pt_end)];
    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(proxy_layer_size^2, ...
        x_pt_start - proxy_layer_size * hx, x_pt_start - hx, ...
        y_pt_start - proxy_layer_size * hy, y_pt_start - hy)];
    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(proxy_layer_size^2, ...
        x_pt_start - proxy_layer_size * hx, x_pt_start - hx, ...
        y_pt_end + hy, y_pt_end + proxy_layer_size * hy)];
    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(proxy_layer_size^2, ...
        x_pt_end + hx, x_pt_end + proxy_layer_size * hx, ...
        y_pt_start - proxy_layer_size * hy, y_pt_start - hy)];
    col_proxy_surface_real = [col_proxy_surface_real; ...
        RandRectangular(proxy_layer_size^2, ...
        x_pt_end + hx, x_pt_end + proxy_layer_size * hx, ...
        y_pt_end + hy, y_pt_end + proxy_layer_size * hy)];

    tmp_ind = col_proxy_surface_real(:, 1) < 0;
    col_proxy_surface_real(tmp_ind, 1) = col_proxy_surface_real(tmp_ind, 1) + 1;
    tmp_ind = col_proxy_surface_real(:, 1) >= 1;
    col_proxy_surface_real(tmp_ind, 1) = col_proxy_surface_real(tmp_ind, 1) - 1;

    tmp_ind = col_proxy_surface_real(:, 2) < 0;
    col_proxy_surface_real(tmp_ind, 2) = col_proxy_surface_real(tmp_ind, 2) + 1;
    tmp_ind = col_proxy_surface_real(:, 2) >= 1;
    col_proxy_surface_real(tmp_ind, 2) = col_proxy_surface_real(tmp_ind, 2) - 1;
    
    col_proxy_surface = [...
        exp(-2 * pi * 1i * col_proxy_surface_real(:, 1)), ...
        exp(-2 * pi * 1i * col_proxy_surface_real(:, 2))];

    % figure();
    % plot(A.row_xy_(:, 1), A.row_xy_(:, 2), "rx", "DisplayName", "row pts");
    % hold on;
    % plot(row_proxy_surface_real(:, 1), row_proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
    % legend;
    % xlim([0, 1]);
    % ylim([0, 1]);
    % axis equal;

    % figure();
    % plot(A.col_pos_(:, 1) / nx, A.col_pos_(:, 2) / ny, "rx", "DisplayName", "col pts");
    % hold on;
    % plot(col_proxy_surface_real(:, 1), col_proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
    % legend;
    % xlim([0, 1]);
    % ylim([0, 1]);
    % axis equal;

    % Construct U using proxy surface.
    gamma_x_I = exp(-2 * pi * 1i * A.row_xy_(:, 1));
    gamma_y_I = exp(-2 * pi * 1i * A.row_xy_(:, 2));
    Ax_I_proxy = x_kernel_fun(gamma_x_I, row_proxy_surface(:, 1));
    Ay_I_proxy = y_kernel_fun(gamma_y_I, row_proxy_surface(:, 2));
    A_I_proxy = Ax_I_proxy .* Ay_I_proxy;
    [row_sk, U, A.row_rank_, ~] = LowRank_Row_ID(A_I_proxy, rank_or_tol);

    % Construct V using proxy surface.
    xi_x_J = exp(-2 * pi * 1i * A.col_pos_(:, 1) / nx);
    xi_y_J = exp(-2 * pi * 1i * A.col_pos_(:, 2) / ny);
    Ax_proxy_J = x_kernel_fun(col_proxy_surface(:, 1), xi_x_J);
    Ay_proxy_J = y_kernel_fun(col_proxy_surface(:, 2), xi_y_J);
    A_proxy_J = Ax_proxy_J .* Ay_proxy_J;
    [col_sk, V, A.col_rank_, ~] = LowRank_ID(A_proxy_J, rank_or_tol);

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