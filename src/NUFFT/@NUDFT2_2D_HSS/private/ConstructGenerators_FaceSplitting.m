function data = ConstructGenerators_FaceSplitting(A, Ax, Ay, Fx, Fy, index, is_root)
% Recursively assemble raw generators on the prebuilt product tree.

nx = Ax.col_size_;
ny = Ay.col_size_;
rx = Ax.row_rank_ * (Ax.level_ > 0);
ry = Ay.row_rank_ * (Ay.level_ > 0);
vx = Ax.col_rank_ * (Ax.level_ > 0);
vy = Ay.col_rank_ * (Ay.level_ > 0);
data.nx = nx;
data.ny = ny;
% Boundaries of the y-only, x-only and both-axis interaction groups.
data.row_ranges = [0, nx * ry, nx * ry + rx * ny, nx * ry + rx * ny + rx * ry];
data.col_ranges = [0, nx * vy, nx * vy + vx * ny, nx * vy + vx * ny + vx * vy];

if A.leaf_ == 1
    p_leaf = index.p(A.row_offset_ + (1 : A.row_size_));
    row_x = index.posx(p_leaf) - Ax.row_offset_;
    row_y = index.posy(p_leaf) - Ay.row_offset_;
    M = A.row_size_;
    Dx = reshape(Ax.Amat_(row_x, :), M, nx, 1);
    Dy = reshape(Ay.Amat_(row_y, :), M, 1, ny);
    A.Amat_ = reshape(Dx .* Dy, M, nx * ny);
    if is_root
        return;
    end

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
    Ux = reshape(Ux, M, rx, 1);
    Uy = reshape(Uy, M, 1, ry);
    U = zeros(M, data.row_ranges(4));
    V = zeros(nx * ny, data.col_ranges(4));
    U(:, 1 : data.row_ranges(2)) = reshape(Dx .* Uy, M, nx * ry);
    U(:, (data.row_ranges(2) + 1) : data.row_ranges(3)) = reshape(Ux .* Dy, M, rx * ny);
    U(:, (data.row_ranges(3) + 1) : data.row_ranges(4)) = reshape(Ux .* Uy, M, rx * ry);
    V(:, 1 : data.col_ranges(2)) = kron(Vy, eye(nx));
    V(:, (data.col_ranges(2) + 1) : data.col_ranges(3)) = kron(eye(ny), Vx);
    V(:, (data.col_ranges(3) + 1) : data.col_ranges(4)) = kron(Vy, Vx);
else
    num_x = max(Ax.num_children_, 1);
    children = cell(1, A.num_children_);
    row_size = 0;
    col_size = 0;
    for i = 1 : A.num_children_
        ix = mod(i - 1, num_x) + 1;
        iy = floor((i - 1) / num_x) + 1;
        if Ax.leaf_ == 1
            children{i} = ConstructGenerators_FaceSplitting( ...
                A.children_{i}, Ax, Ay.children_{iy}, Fx, Fy.children{iy}, index, false);
        elseif Ay.leaf_ == 1
            children{i} = ConstructGenerators_FaceSplitting( ...
                A.children_{i}, Ax.children_{ix}, Ay, Fx.children{ix}, Fy, index, false);
        else
            children{i} = ConstructGenerators_FaceSplitting( ...
                A.children_{i}, Ax.children_{ix}, Ay.children_{iy}, ...
                Fx.children{ix}, Fy.children{iy}, index, false);
        end
        row_size = row_size + A.children_{i}.row_rank_;
        col_size = col_size + A.children_{i}.col_rank_;
    end

    A.Bmat_ = cell(A.num_children_);
    for i = 1 : A.num_children_
        ix = mod(i - 1, num_x) + 1;
        iy = floor((i - 1) / num_x) + 1;
        for j = 1 : A.num_children_
            if i == j
                continue;
            end
            jx = mod(j - 1, num_x) + 1;
            jy = floor((j - 1) / num_x) + 1;
            if ix == jx
                g = 1;
                K = kron(Ay.Bmat_{iy, jy}, eye(children{i}.nx));
            elseif iy == jy
                g = 2;
                K = kron(eye(children{i}.ny), Ax.Bmat_{ix, jx});
            else
                g = 3;
                K = kron(Ay.Bmat_{iy, jy}, Ax.Bmat_{ix, jx});
            end
            I = (children{i}.row_ranges(g) + 1) : children{i}.row_ranges(g + 1);
            J = (children{j}.col_ranges(g) + 1) : children{j}.col_ranges(g + 1);
            A.Bmat_{i, j} = zeros(A.children_{i}.row_rank_, A.children_{j}.col_rank_);
            A.Bmat_{i, j}(I, J) = K;
        end
    end
    if is_root
        return;
    end

    U = zeros(row_size, data.row_ranges(4));
    V = zeros(col_size, data.col_ranges(4));
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        ix = mod(i - 1, num_x) + 1;
        iy = floor((i - 1) / num_x) + 1;
        if Ax.leaf_ == 1
            Ex = eye(nx); Rx = eye(rx); Wx = eye(vx);
            Fxi = zeros(rx, nx);
        else
            Ex = zeros(Ax.children_{ix}.col_size_, nx);
            J = Ax.children_{ix}.col_offset_ - Ax.col_offset_ + (1 : Ax.children_{ix}.col_size_);
            Ex(:, J) = eye(Ax.children_{ix}.col_size_);
            Rx = Ax.children_{ix}.Rmat_; Wx = Ax.children_{ix}.Wmat_; Fxi = Fx.mat{ix};
        end
        if Ay.leaf_ == 1
            Ey = eye(ny); Ry = eye(ry); Wy = eye(vy);
            Fyi = zeros(ry, ny);
        else
            Ey = zeros(Ay.children_{iy}.col_size_, ny);
            J = Ay.children_{iy}.col_offset_ - Ay.col_offset_ + (1 : Ay.children_{iy}.col_size_);
            Ey(:, J) = eye(Ay.children_{iy}.col_size_);
            Ry = Ay.children_{iy}.Rmat_; Wy = Ay.children_{iy}.Wmat_; Fyi = Fy.mat{iy};
        end
        % Fill the nonzero blocks of the three interaction groups.
        I = row_offset + (1 : children{i}.row_ranges(2));
        U(I, 1 : data.row_ranges(2)) = kron(Ry, Ex);
        I = row_offset + ((children{i}.row_ranges(2) + 1) : children{i}.row_ranges(3));
        U(I, (data.row_ranges(2) + 1) : data.row_ranges(3)) = kron(Ey, Rx);
        I = row_offset + ((children{i}.row_ranges(3) + 1) : children{i}.row_ranges(4));
        U(I, 1 : data.row_ranges(2)) = kron(Ry, Fxi);
        U(I, (data.row_ranges(2) + 1) : data.row_ranges(3)) = kron(Fyi, Rx);
        U(I, (data.row_ranges(3) + 1) : data.row_ranges(4)) = kron(Ry, Rx);
        J = col_offset + (1 : children{i}.col_ranges(2));
        V(J, 1 : data.col_ranges(2)) = kron(Wy, Ex);
        J = col_offset + ((children{i}.col_ranges(2) + 1) : children{i}.col_ranges(3));
        V(J, (data.col_ranges(2) + 1) : data.col_ranges(3)) = kron(Ey, Wx);
        J = col_offset + ((children{i}.col_ranges(3) + 1) : children{i}.col_ranges(4));
        V(J, (data.col_ranges(3) + 1) : data.col_ranges(4)) = kron(Wy, Wx);
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end
end

A.row_rank_ = size(U, 2);
A.col_rank_ = size(V, 2);
if A.leaf_ == 1
    A.Umat_ = U;
    A.Vmat_ = V;
else
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        A.children_{i}.Rmat_ = U(row_offset + (1 : A.children_{i}.row_rank_), :);
        A.children_{i}.Wmat_ = V(col_offset + (1 : A.children_{i}.col_rank_), :);
        row_offset = row_offset + A.children_{i}.row_rank_;
        col_offset = col_offset + A.children_{i}.col_rank_;
    end
end

end
