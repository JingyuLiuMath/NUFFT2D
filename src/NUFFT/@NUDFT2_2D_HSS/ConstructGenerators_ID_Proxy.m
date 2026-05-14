function ConstructGenerators_ID_Proxy(A, level, tol)

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
    tol (1, 1) double;
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

    % if A.leaf_ == 1
    %     tmp_ind = proxy_surface_real(:, 1) < 0;
    %     proxy_surface_real(tmp_ind, 1) = proxy_surface_real(tmp_ind, 1) + 1;
    %     tmp_ind = proxy_surface_real(:, 1) >= 1;
    %     proxy_surface_real(tmp_ind, 1) = proxy_surface_real(tmp_ind, 1) - 1;
    % 
    %     tmp_ind = proxy_surface_real(:, 2) < 0;
    %     proxy_surface_real(tmp_ind, 2) = proxy_surface_real(tmp_ind, 2) + 1;
    %     tmp_ind = proxy_surface_real(:, 2) >= 1;
    %     proxy_surface_real(tmp_ind, 2) = proxy_surface_real(tmp_ind, 2) - 1;
    % 
    % 
    %     figure();
    %     plot(A.row_xy_(:, 1), A.row_xy_(:, 2), "rx", "DisplayName", "row pts");
    %     hold on;
    %     plot(proxy_surface_real(:, 1), proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
    %     legend;
    %     xlim([0, 1]);
    %     ylim([0, 1]);
    %     axis equal;
    %     filename = "./figure/row_proxy_surface_" ...
    %         + string(A.level_) ...
    %         + "_" + string(A.row_offset_)  ...
    %         + "_" + string(A.col_offset_)  ...
    %         + "_" + string(A.row_size_) ...
    %         + "_" + string(A.col_size_);
    %     saveas(gcf, filename + ".png", "png");
    %     saveas(gcf, filename + ".eps", "epsc");
    % 
    %     figure();
    %     plot(A.col_pos_(:, 1) / nx, A.col_pos_(:, 2) / ny, "rx", "DisplayName", "col pts");
    %     hold on;
    %     plot(proxy_surface_real(:, 1), proxy_surface_real(:, 2), "bx", "DisplayName", "proxy pts");
    %     legend;
    %     xlim([0, 1]);
    %     ylim([0, 1]);
    %     axis equal;
    %     filename = "./figure/col_proxy_surface_" ...
    %         + string(A.level_) ...
    %         + "_" + string(A.row_offset_)  ...
    %         + "_" + string(A.col_offset_)  ...
    %         + "_" + string(A.row_size_) ...
    %         + "_" + string(A.col_size_);
    %     saveas(gcf, filename + ".png", "png");
    %     saveas(gcf, filename + ".eps", "epsc");
    % end

    % Construct U using proxy surface.
    gamma_x_I = exp(-2 * pi * 1i * A.row_xy_(:, 1));
    gamma_y_I = exp(-2 * pi * 1i * A.row_xy_(:, 2));
    Ax_I_proxy = x_kernel_fun(gamma_x_I, proxy_surface(:, 1));
    Ay_I_proxy = y_kernel_fun(gamma_y_I, proxy_surface(:, 2));
    A_I_proxy = Ax_I_proxy .* Ay_I_proxy;
    [row_sk, U, A.row_rank_, ~] = LowRank_Row_ID(A_I_proxy, tol);

    % Construct V using proxy surface.
    xi_x_J = exp(-2 * pi * 1i * A.col_pos_(:, 1) / nx);
    xi_y_J = exp(-2 * pi * 1i * A.col_pos_(:, 2) / ny);
    Ax_proxy_J = x_kernel_fun(proxy_surface(:, 1), xi_x_J);
    Ay_proxy_J = y_kernel_fun(proxy_surface(:, 2), xi_y_J);
    A_proxy_J = Ax_proxy_J .* Ay_proxy_J;
    [col_sk, V, A.col_rank_, ~] = LowRank_ID(A_proxy_J, tol);

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
        A.children_{i}.ConstructGenerators_ID_Proxy(level, tol);
    end
end

end