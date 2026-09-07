function ConstructGenerators_FaceSplitting(A, level, Ax, Ay, px, py, Fx, Fy)

arguments (Input)
    A NUDFT2_2D_HSS;
    level (1, 1) double;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    px (:, 1) double;
    py (:, 1) double;
    Fx struct;
    Fy struct;
end

if A.leaf_ == 0
    xc = Ax.children_;
    yc = Ay.children_;
    fxc = Fx.children;
    fyc = Fy.children;
    if Ax.leaf_
        xc = {Ax};
        fxc = {Fx};
    end
    if Ay.leaf_
        yc = {Ay};
        fyc = {Fy};
    end
    [ix, iy] = ndgrid(1 : numel(xc), 1 : numel(yc));
    ix = ix(:);
    iy = iy(:);
end

if A.level_ == level
    nx = Ax.col_size_;
    ny = Ay.col_size_;
    % A source root has no basis for interactions outside the source.
    row_rank_x = Ax.row_rank_ * (Ax.level_ > 0);
    row_rank_y = Ay.row_rank_ * (Ay.level_ > 0);
    col_rank_x = Ax.col_rank_ * (Ax.level_ > 0);
    col_rank_y = Ay.col_rank_ * (Ay.level_ > 0);
    if A.leaf_ == 1
        row_x = px(Ax.row_offset_ + (1 : Ax.row_size_));
        row_y = py(Ay.row_offset_ + (1 : Ay.row_size_));
        [~, row_x, row_y] = intersect(row_x, row_y, "sorted");
        M = A.row_size_;
        Dx = reshape(Ax.Amat_(row_x, :), M, nx, 1);
        Dy = reshape(Ay.Amat_(row_y, :), M, 1, ny);
        A.Amat_ = reshape(Dx .* Dy, M, nx * ny);

        if Ax.level_ == 0
            Ux = zeros(M, 0);
            Vx = zeros(nx, 0);
        else
            Ux = Ax.Umat_(row_x, :);
            Vx = Ax.Vmat_;
        end
        if Ay.level_ == 0
            Uy = zeros(M, 0);
            Vy = zeros(ny, 0);
        else
            Uy = Ay.Umat_(row_y, :);
            Vy = Ay.Vmat_;
        end
        % The three groups describe y, x, and both-axis interactions.
        Ux = reshape(Ux, M, row_rank_x, 1);
        Uy = reshape(Uy, M, 1, row_rank_y);
        U = [reshape(Dx .* Uy, M, nx * row_rank_y), ...
            reshape(Ux .* Dy, M, row_rank_x * ny), ...
            reshape(Ux .* Uy, M, row_rank_x * row_rank_y)];
        V = [kron(Vy, eye(nx)), kron(eye(ny), Vx), kron(Vy, Vx)];
    else
        row_ranges = zeros(A.num_children_, 4);
        col_ranges = zeros(A.num_children_, 4);
        row_size = 0;
        col_size = 0;
        for i = 1 : A.num_children_
            Cx = xc{ix(i)};
            Cy = yc{iy(i)};
            rx = Cx.row_rank_ * (Cx.level_ > 0);
            ry = Cy.row_rank_ * (Cy.level_ > 0);
            vx = Cx.col_rank_ * (Cx.level_ > 0);
            vy = Cy.col_rank_ * (Cy.level_ > 0);
            row_ranges(i, :) = cumsum([0, Cx.col_size_ * ry, rx * Cy.col_size_, rx * ry]);
            col_ranges(i, :) = cumsum([0, Cx.col_size_ * vy, vx * Cy.col_size_, vx * vy]);
            row_size = row_size + A.children_{i}.row_rank_;
            col_size = col_size + A.children_{i}.col_rank_;
        end
        U = zeros(row_size, nx * row_rank_y + row_rank_x * ny + row_rank_x * row_rank_y);
        V = zeros(col_size, nx * col_rank_y + col_rank_x * ny + col_rank_x * col_rank_y);

        A.Bmat_ = cell(A.num_children_);
        row_offset = 0;
        col_offset = 0;
        for i = 1 : A.num_children_
            for j = [1 : (i - 1), (i + 1) : A.num_children_]
                if ix(i) == ix(j)
                    g = 1;
                    K = kron(Ay.Bmat_{iy(i), iy(j)}, eye(xc{ix(i)}.col_size_));
                elseif iy(i) == iy(j)
                    g = 2;
                    K = kron(eye(yc{iy(i)}.col_size_), Ax.Bmat_{ix(i), ix(j)});
                else
                    g = 3;
                    K = kron(Ay.Bmat_{iy(i), iy(j)}, Ax.Bmat_{ix(i), ix(j)});
                end
                I = (row_ranges(i, g) + 1) : row_ranges(i, g + 1);
                J = (col_ranges(j, g) + 1) : col_ranges(j, g + 1);
                A.Bmat_{i, j} = A.children_{i}.FaceSplitting_Tu_(:, I) ...
                    * K * A.children_{j}.FaceSplitting_Tv_(:, J)';
            end
            % E follows the column partition; R and W are source generators.
            Cx = xc{ix(i)};
            Cy = yc{iy(i)};
            if Ax.leaf_
                Ex = eye(nx);
                Rx = eye(row_rank_x);
                Wx = eye(col_rank_x);
                Fxi = zeros(row_rank_x, nx);
            else
                Ex = zeros(Cx.col_size_, nx);
                J = Cx.col_offset_ - Ax.col_offset_ + (1 : Cx.col_size_);
                Ex(:, J) = eye(Cx.col_size_);
                Rx = Ax.Rmat_{ix(i)};
                Wx = Ax.Wmat_{ix(i)};
                Fxi = Fx.mat{ix(i)};
            end
            if Ay.leaf_
                Ey = eye(ny);
                Ry = eye(row_rank_y);
                Wy = eye(col_rank_y);
                Fyi = zeros(row_rank_y, ny);
            else
                Ey = zeros(Cy.col_size_, ny);
                J = Cy.col_offset_ - Ay.col_offset_ + (1 : Cy.col_size_);
                Ey(:, J) = eye(Cy.col_size_);
                Ry = Ay.Rmat_{iy(i)};
                Wy = Ay.Wmat_{iy(i)};
                Fyi = Fy.mat{iy(i)};
            end
            C = A.children_{i};
            Tu1 = C.FaceSplitting_Tu_(:, 1 : row_ranges(i, 2));
            Tu2 = C.FaceSplitting_Tu_(:, (row_ranges(i, 2) + 1) : row_ranges(i, 3));
            Tu3 = C.FaceSplitting_Tu_(:, (row_ranges(i, 3) + 1) : row_ranges(i, 4));
            Tv1 = C.FaceSplitting_Tv_(:, 1 : col_ranges(i, 2));
            Tv2 = C.FaceSplitting_Tv_(:, (col_ranges(i, 2) + 1) : col_ranges(i, 3));
            Tv3 = C.FaceSplitting_Tv_(:, (col_ranges(i, 3) + 1) : col_ranges(i, 4));
            % Multiply the three transfer blocks directly.
            U(row_offset + (1 : C.row_rank_), :) = [ ...
                Tu1 * kron(Ry, Ex) + Tu3 * kron(Ry, Fxi), ...
                Tu2 * kron(Ey, Rx) + Tu3 * kron(Fyi, Rx), ...
                Tu3 * kron(Ry, Rx)];
            V(col_offset + (1 : C.col_rank_), :) = [ ...
                Tv1 * kron(Wy, Ex), ...
                Tv2 * kron(Ey, Wx), ...
                Tv3 * kron(Wy, Wx)];
            row_offset = row_offset + C.row_rank_;
            col_offset = col_offset + C.col_rank_;
        end
        % Clear after all couplings and transfers have used the child factors.
        for i = 1 : A.num_children_
            A.children_{i}.FaceSplitting_Tu_ = [];
            A.children_{i}.FaceSplitting_Tv_ = [];
        end
    end

    [U, A.FaceSplitting_Tu_] = qr(U, "econ");
    [V, A.FaceSplitting_Tv_] = qr(V, "econ");
    A.row_rank_ = size(U, 2);
    A.col_rank_ = size(V, 2);
    if A.leaf_ == 1
        A.Umat_ = U;
        A.Vmat_ = V;
    else
        row_offset = 0;
        col_offset = 0;
        for i = 1 : A.num_children_
            C = A.children_{i};
            A.Rmat_{i} = U(row_offset + (1 : C.row_rank_), :);
            A.Wmat_{i} = V(col_offset + (1 : C.col_rank_), :);
            row_offset = row_offset + C.row_rank_;
            col_offset = col_offset + C.col_rank_;
        end
    end
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.ConstructGenerators_FaceSplitting( ...
            level, xc{ix(i)}, yc{iy(i)}, px, py, fxc{ix(i)}, fyc{iy(i)});
    end
end

end
