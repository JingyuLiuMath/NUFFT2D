function ConstructGenerators_ID_Proxy0(A, level, tol, ...
    xy, col_pos, ...
    row_inact, col_inact)

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
    tol (1, 1) double;
    xy (:, 2) double;
    col_pos (:, 2) double;
    row_inact (:, 1) double;
    col_inact (:, 1) double;
end

M = A.row_global_size_;
N = A.col_global_size_;

nx = A.nx_;
ny = A.ny_;
hx = 1 / nx;
hhx = hx / 2;
hy = 1 / ny;
hhy = hy / 2;
n = max(nx, ny);

x_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, nx);
y_kernel_fun = @(z, w) NUFFT2_Kernel(z, w, ny);

local_nx = A.x_col_size_;
local_ny = A.y_col_size_;
local_n = max(local_nx, local_ny);

proxy_layer_size = 2 + log2(local_n);
rank_1d = ceil(log(4 / tol) * log(2 * pi * (2 * local_n + 1)) * 2 / pi^2);
sampling_size_cross = ceil(1.2 * rank_1d * local_n);
sampling_size_diag = rank_1d^2;

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

    % Generate proxy surface.
    proxy_surface_real = [];

    x_pt_start = A.x_pos_start_ / nx;
    x_pt_end = A.x_pos_end_ / nx ;
    y_pt_start = A.y_pos_start_ / ny;
    y_pt_end = A.y_pos_end_ / ny;

    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_cross, ...
        x_pt_start - hhx, x_pt_end + hhx, ...
        y_pt_start - hhy - proxy_layer_size * hy, y_pt_start - hhy)];
    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_cross, ...
        x_pt_start - hhx, x_pt_end + hhx, ...
        y_pt_end + hhy, y_pt_end + hhy + proxy_layer_size * hy)];
    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_cross, ...
        x_pt_start - hhx - proxy_layer_size * hx, x_pt_start - hhx, ...
        y_pt_start - hhy, y_pt_end + hhy)];
    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_cross, ...
        x_pt_end + hhx, x_pt_end + hhx + proxy_layer_size * hx, ...
        y_pt_start - hhy, y_pt_end + hhy)];

    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_diag, ...
        x_pt_start - hhx - proxy_layer_size * hx, x_pt_start - hhx, ...
        y_pt_start - hhy - proxy_layer_size * hy, y_pt_start - hhy)];
    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_diag, ...
        x_pt_start - hhx - proxy_layer_size * hx, x_pt_start - hhx, ...
        y_pt_end + hhy, y_pt_end + hhy + proxy_layer_size * hy)];
    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_diag, ...
        x_pt_end + hhx, x_pt_end + hhx + proxy_layer_size * hx, ...
        y_pt_start - hhy - proxy_layer_size * hy, y_pt_start - hhy)];
    proxy_surface_real = [proxy_surface_real; ...
        RandRectangular(sampling_size_diag, ...
        x_pt_end + hhx, x_pt_end + hhx + proxy_layer_size * hx, ...
        y_pt_end + hhy, y_pt_end + hhy + proxy_layer_size * hy)];


    proxy_surface = [...
        exp(-2 * pi * 1i * proxy_surface_real(:, 1)), ...
        exp(-2 * pi * 1i * proxy_surface_real(:, 2))];

    % Construct U using proxy surface.
    xy_I = xy(A.row_ind_, :);
    gamma_x_I = exp(-2 * pi * 1i * xy_I(:, 1));
    gamma_y_I = exp(-2 * pi * 1i * xy_I(:, 2));
    Ax_I_proxy = x_kernel_fun(gamma_x_I, proxy_surface(:, 1));
    Ay_I_proxy = y_kernel_fun(gamma_y_I, proxy_surface(:, 2));
    A_I_proxy = Ax_I_proxy .* Ay_I_proxy;
    [row_sk, U, A.row_rank_, row_re] = LowRank_Row_ID(A_I_proxy, tol);

    E_I_proxy = A_I_proxy - U * A_I_proxy(row_sk, :);
    rel_err_I_proxy = norm(E_I_proxy, "fro") / norm(A_I_proxy, "fro");
    fprintf("rel_err_I_proxy: %.4e\n", rel_err_I_proxy);

    Jc = (1 : N)';
    self_col_ind = A.col_offset_ + (1 : A.col_size_)';
    Jc([self_col_ind; col_inact]) = [];
    col_pos_Jc = col_pos(Jc, :);
    xi_x_Jc = exp(-2 * pi * 1i * col_pos_Jc(:, 1) / nx);
    xi_y_Jc = exp(-2 * pi * 1i * col_pos_Jc(:, 2) / ny);
    Ax_I_Jc = x_kernel_fun(gamma_x_I, xi_x_Jc);
    Ay_I_Jc = y_kernel_fun(gamma_y_I, xi_y_Jc);
    A_I_Jc = Ax_I_Jc .* Ay_I_Jc;
    
    E_I_Jc = A_I_Jc - U * A_I_Jc(row_sk, :);
    rel_err_I_Jc = norm(E_I_Jc, "fro") / norm(A_I_Jc, "fro");
    fprintf("rel_err_I_Jc: %.4e\n", rel_err_I_Jc);

    if rel_err_I_Jc > tol * 10
        figure();
        plot(xy_I(:, 1), xy_I(:, 2), "rx", "DisplayName", "row pts");
        hold on;
        plot(col_pos_Jc(:, 1) / nx, col_pos_Jc(:, 2) / ny, "go", "DisplayName", "HSS row pts");
        plot(proxy_surface_real(:, 1), proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
        legend;
        xlim([0, 1]);
        ylim([0, 1]);
        axis equal;

        [tmp_row_sk, tmp_U, tmp_row_rank, tmp_row_re] = LowRank_Row_ID(A_I_Jc, tol);

        tmp_E_I_Jc = A_I_Jc - tmp_U * A_I_Jc(tmp_row_sk, :);
        tmp_rel_err_I_Jc = norm(tmp_E_I_Jc, "fro") / norm(A_I_Jc, "fro");
        fprintf("tmp_rel_err_I_Jc: %.4e\n", tmp_rel_err_I_Jc);

        tmp_E_I_proxy = A_I_proxy - tmp_U * A_I_proxy(tmp_row_sk, :);
        tmp_rel_err_I_proxy = norm(tmp_E_I_proxy, "fro") / norm(A_I_proxy, "fro");
        fprintf("tmp_rel_err_I_proxy: %.4e\n", tmp_rel_err_I_proxy);

        [tmp_col_sk, tmp_V, ~, ~] = LowRank_ID(A_I_Jc, tol);
        figure();
        plot(xy_I(:, 1), xy_I(:, 2), "rx", "DisplayName", "row pts");
        hold on;
        plot(col_pos_Jc(tmp_col_sk, 1) / nx, col_pos_Jc(tmp_col_sk, 2) / ny, "go", "DisplayName", "HSS col pts");
        plot(proxy_surface_real(:, 1), proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
        legend;
        xlim([0, 1]);
        ylim([0, 1]);
        axis equal;

        keyboard;
    end

    % Construct V using proxy surface.
    col_pos_J = col_pos(A.col_ind_, :);
    xi_x_J = exp(-2 * pi * 1i * col_pos_J(:, 1) / nx);
    xi_y_J = exp(-2 * pi * 1i * col_pos_J(:, 2) / ny);
    Ax_proxy_J = x_kernel_fun(proxy_surface(:, 1), xi_x_J);
    Ay_proxy_J = y_kernel_fun(proxy_surface(:, 2), xi_y_J);
    A_proxy_J = Ax_proxy_J .* Ay_proxy_J;
    [col_sk, V, A.col_rank_, col_re] = LowRank_ID(A_proxy_J, tol);
    
    % figure();
    % plot(col_pos_J(:, 1) / nx, col_pos_J(:, 2) / ny, "go", "DisplayName", "col pts");
    % hold on;
    % plot(col_proxy_surface_real(:, 1), col_proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
    % legend;
    % xlim([0, 1]);
    % ylim([0, 1]);
    % axis equal;

    E_proxy_J = A_proxy_J - A_proxy_J(:, col_sk) * V';
    rel_err_proxy_J = norm(E_proxy_J, "fro") / norm(A_proxy_J, "fro");
    fprintf("rel_err_proxy_J: %.4e\n", rel_err_proxy_J);

    Ic = (1 : M)';
    self_row_ind = A.row_offset_ + (1 : A.row_size_)';
    Ic([self_row_ind; row_inact]) = [];
    xy_Ic = xy(Ic, :);
    gamma_x_Ic = exp(-2 * pi * 1i * xy_Ic(:, 1));
    gamma_y_Ic = exp(-2 * pi * 1i * xy_Ic(:, 2));
    Ax_Ic_J = x_kernel_fun(gamma_x_Ic, xi_x_J);
    Ay_Ic_J = y_kernel_fun(gamma_y_Ic, xi_y_J);
    A_Ic_J = Ax_Ic_J .* Ay_Ic_J;

    E_Ic_J = A_Ic_J - A_Ic_J(:, col_sk) * V';
    rel_err_Ic_J = norm(E_Ic_J, "fro") / norm(A_Ic_J, "fro");
    fprintf("rel_err_Ic_J: %.4e\n", rel_err_Ic_J);

    if rel_err_Ic_J > tol * 10
        figure();
        plot(col_pos_J(:, 1) / nx, col_pos_J(:, 2) / ny, "rx", "DisplayName", "col pts");
        hold on;
        plot(xy_Ic(:, 1), xy_Ic(:, 2), "go", "DisplayName", "HSS col pts");
        plot(proxy_surface_real(:, 1), proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
        legend;
        xlim([0, 1]);
        ylim([0, 1]);
        axis equal;

        [tmp_col_sk, tmp_V, tmp_col_rank, tmp_col_re] = LowRank_ID(A_Ic_J, tol);

        tmp_E_proxy_J = A_proxy_J - A_proxy_J(:, tmp_col_sk) * tmp_V';
        tmp_rel_err_proxy_J = norm(tmp_E_proxy_J, "fro") / norm(A_proxy_J, "fro");
        fprintf("tmp_rel_err_proxy_J: %.4e\n", tmp_rel_err_proxy_J);

        tmp_E_Ic_J = A_Ic_J - A_Ic_J(:, tmp_col_sk) * tmp_V';
        tmp_rel_err_Ic_J = norm(tmp_E_Ic_J, "fro") / norm(A_Ic_J, "fro");
        fprintf("tmp_rel_err_Ic_J: %.4e\n", tmp_rel_err_Ic_J);

        [tmp_row_sk, tmp_U, ~, ~] = LowRank_ID(A_Ic_J, tol);
        figure();
        plot(col_pos_J(:, 1) / nx, col_pos_J(:, 2) / ny, "rx", "DisplayName", "col pts");
        hold on;
        plot(xy_Ic(tmp_row_sk, 1), xy_Ic(tmp_row_sk, 2), "go", "DisplayName", "HSS row pts");
        plot(proxy_surface_real(:, 1), proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
        legend;
        xlim([0, 1]);
        ylim([0, 1]);
        axis equal;

        keyboard;
    end

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
            A.children_{i}.row_xy_ = [];
            A.children_{i}.col_pos_ = [];
        end
    end

    % Update row and col.
    A.row_re_ = A.row_ind_(row_re);
    A.col_re_ = A.col_ind_(col_re);
    A.row_ind_ = A.row_ind_(row_sk);
    A.col_ind_ = A.col_ind_(col_sk);
else
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_ID_Proxy0(level, tol, ...
            xy, col_pos, row_inact, col_inact);
    end
end

end