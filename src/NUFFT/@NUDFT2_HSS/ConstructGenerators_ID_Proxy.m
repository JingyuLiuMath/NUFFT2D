function ConstructGenerators_ID_Proxy(A, level, rank_or_tol)
% ConstructGenerators

% Jingyu Liu, November 24, 2024.

arguments (Input)
    A NUDFT2_HSS;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

N = A.col_global_size_;
half_length = 0.5 / N;
kernel_fun = @(z, w) NUFFT2_Kernel(z, w, N);

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

    row_x = A.row_x_;
    tmp_ind = find(row_x >= (1 - 1 / N / 2));
    row_x(tmp_ind) = row_x(tmp_ind) - 1;
    if rank_or_tol >= 1
        rank_or_tol = min([A.col_size_, rank_or_tol]);
        proxy_size = min([A.col_size_, 2 * rank_or_tol]);
    else
        Ir = [...
            exp(-2 * pi * 1i * min(row_x)), ...
            exp(-2 * pi * 1i * max(row_x)), ...
            exp(-2 * pi * 1i * (A.pos_end_ + 1) / N), ...
            exp(-2 * pi * 1i * (A.pos_start_ - 1) / N)];
        [~, ~, ~, cr] = mobiusT(Ir);
        k = ceil(1/pi^2*log(4/rank_or_tol)*log(16*cr));

        Ic = [...
            exp(-2 * pi * 1i * (A.pos_end_ + 1 / 2) / N), ...
            exp(-2 * pi * 1i * (A.pos_start_ - 1 / 2) / N), ...
            exp(-2 * pi * 1i * min(A.col_pos_) / N), ...
            exp(-2 * pi * 1i * max(A.col_pos_) / N)];
        [~, ~, ~, cc] = mobiusT(Ic);
        s = ceil(1/pi^2*log(4/rank_or_tol)*log(16*cc));
        rank_or_tol = max(k, s);

        proxy_size = min([A.col_size_, rank_or_tol]);
    end

    % Generate proxy surface.
    col_x_start = A.pos_start_ / N;
    col_x_end = A.pos_end_ / N;
    x_prony = GenerateProxySurface(N, ...
        col_x_start - half_length, ...
        col_x_end + half_length, ...
        proxy_size);
    proxy_surface = exp(-2 * pi * 1i * x_prony);

    % Construct U using proxy surface.
    z_I = exp(-2 * pi * 1i * A.row_x_);  % row points.
    A_I_proxy = kernel_fun(z_I, proxy_surface);
    [row_sk, U, A.row_rank_] ...
        = LowRank_Row_ID(A_I_proxy, rank_or_tol);
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
    w_J = exp(-2 * pi * 1i * A.col_pos_ / N);  % col points.
    A_proxy_J = kernel_fun(proxy_surface, w_J);
    [col_sk, V, A.col_rank_] ...
        = LowRank_ID(A_proxy_J, rank_or_tol);

    if A.leaf_ == 1
        % Assign U and V.
        A.Umat_ = U;
        A.Vmat_ = V;

        % Assign full mat.
        A.Amat_ = kernel_fun(z_I, w_J);
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
        A.children_{i}.ConstructGenerators_ID_Proxy(level, rank_or_tol);
    end
end

end