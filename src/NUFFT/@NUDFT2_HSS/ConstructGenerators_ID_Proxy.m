function ConstructGenerators_ID_Proxy(A, level, tol)

arguments (Input)
    A NUDFT2_HSS;
    level (1, 1) double;
    tol (1, 1) double;
end

M = A.row_global_size_;
N = A.col_global_size_;

h = 1 / N;
hh = h / 2;

kernel_fun = @(z, w) NUFFT2_Kernel(z, w, N);

local_n = A.col_size_;

proxy_layer_size = 2 + log2(local_n);
rank_1d = ceil(log(4 / tol) * log(2 * pi * (2 * local_n + 1)) * 2 / pi^2);
sampling_size = ceil(1.2 * rank_1d);

if A.level_ == level
    if A.leaf_ == 0
        % Merge data from children and construct B.
        A.row_x_ = [];
        A.col_pos_ = [];
        for i = 1 : A.num_children_
            A.row_x_ = [A.row_x_; A.children_{i}.row_x_];
            A.col_pos_ = [A.col_pos_; A.children_{i}.col_pos_];
            c_z_I = exp(-2 * pi * 1i * A.children_{i}.row_x_);
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                c_w_J = exp(-2 * pi * 1i * A.children_{j}.col_pos_ / N);
                A.Bmat_{i, j} = kernel_fun(c_z_I, c_w_J);
            end
        end
    end



    % Generate proxy surface.
    proxy_surface_real = [];

    pt_start = A.pos_start_ / N;
    pt_end = A.pos_end_ / N;

    proxy_surface_real = [proxy_surface_real; ...
        RandInterval(sampling_size, ...
        pt_start - hh - proxy_layer_size * h, pt_start - hh)];
    proxy_surface_real = [proxy_surface_real; ...
        RandInterval(sampling_size, ...
        pt_end + hh, pt_end + hh + proxy_layer_size * h)];

    proxy_surface = exp(-2 * pi * 1i * proxy_surface_real);

    % Construct U using proxy surface.
    gamma_I = exp(-2 * pi * 1i * A.row_x_);  % row points.
    A_I_proxy = kernel_fun(gamma_I, proxy_surface);
    [row_sk, U, A.row_rank_, ~] ...
        = LowRank_Row_ID(A_I_proxy, tol);
    % *****************************************************************
    % Plot row points and proxy surface.
    % figure();
    % finer_N = 4 * N;
    % bg_x = exp(-2 * pi * 1i * (0 : finer_N)' / finer_N);
    % plot(real(bg_x), imag(bg_x), "k-", LineWidth=1);
    % axis equal;
    % hold on;
    % fig_row_points = plot(real(z_I), imag(z_I), ...
    %     Color="#EDB120", Marker="x", ...
    %     LineStyle="none", LineWidth=1);
    % fig_proxy = plot(real(proxy_surface), imag(proxy_surface), ...
    %     Color="#0072BD", Marker="*", ...
    %     LineStyle="none", LineWidth=1);
    % legend([fig_row_points, fig_proxy], ...
    %     "row points", "proxy surface");
    % -----------------------------------------------------------------

    % Construct V using proxy surface.
    xi_J = exp(-2 * pi * 1i * A.col_pos_ / N);  % col points.
    A_proxy_J = kernel_fun(proxy_surface, xi_J);
    [col_sk, V, A.col_rank_, ~] ...
        = LowRank_ID(A_proxy_J, tol);

    if A.leaf_ == 1
        % Assign U and V.
        A.Umat_ = U;
        A.Vmat_ = V;

        % Assign full mat.
        A.Amat_ = kernel_fun(gamma_I, xi_J);
    else
        % Assign R and W.
        offset = 0;
        for i = 1 : A.num_children_
            current_size = size(A.children_{i}.row_x_, 1);
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
            A.children_{i}.row_x_ = [];
            A.children_{i}.col_pos_ = [];
        end
    end

    % Update row and col.
    A.row_x_ = A.row_x_(row_sk);
    A.col_pos_ = A.col_pos_(col_sk);
else
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_ID_Proxy(level, tol);
    end
end

end