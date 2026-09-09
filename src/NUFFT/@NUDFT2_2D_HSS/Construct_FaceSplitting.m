function Construct_FaceSplitting(A, Ax, Ay, Fx, Fy)

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    Fx struct = struct();
    Fy struct = struct();
end

if nargin == 3
    Fx = FaceSplitting_Transfers(Ax);
    Fy = FaceSplitting_Transfers(Ay);
end

num_children_x = max(Ax.num_children_, 1);

% Construct children first.
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    if Ax.leaf_ == 1
        A.children_{i}.Construct_FaceSplitting(...
            Ax, Ay.children_{iy}, Fx, Fy.children{iy});
    elseif Ay.leaf_ == 1
        A.children_{i}.Construct_FaceSplitting(...
            Ax.children_{ix}, Ay, Fx.children{ix}, Fy);
    else
        A.children_{i}.Construct_FaceSplitting(...
            Ax.children_{ix}, Ay.children_{iy}, ...
            Fx.children{ix}, Fy.children{iy});
    end
end

nx = Ax.col_size_;
ny = Ay.col_size_;
row_rank_x = Ax.row_rank_ * (Ax.level_ > 0);
row_rank_y = Ay.row_rank_ * (Ay.level_ > 0);
col_rank_x = Ax.col_rank_ * (Ax.level_ > 0);
col_rank_y = Ay.col_rank_ * (Ay.level_ > 0);
row_ranges = cumsum([0, nx * row_rank_y, ...
    row_rank_x * ny, row_rank_x * row_rank_y]);
col_ranges = cumsum([0, nx * col_rank_y, ...
    col_rank_x * ny, col_rank_x * col_rank_y]);

if A.leaf_ == 1
    M = A.row_size_;
    rows_x = A.row_xy_(:, 1);
    rows_y = A.row_xy_(:, 2);
    Dx = reshape(Ax.Amat_(rows_x, :), M, nx, 1);
    Dy = reshape(Ay.Amat_(rows_y, :), M, 1, ny);
    A.Amat_ = reshape(Dx .* Dy, M, A.col_size_);

    if A.level_ > 0
        if Ax.level_ == 0
            Ux = zeros(M, 0);
            Vx = zeros(nx, 0);
        else
            Ux = Ax.Umat_(rows_x, :);
            Vx = Ax.Vmat_;
        end
        if Ay.level_ == 0
            Uy = zeros(M, 0);
            Vy = zeros(ny, 0);
        else
            Uy = Ay.Umat_(rows_y, :);
            Vy = Ay.Vmat_;
        end

        Ux = reshape(Ux, M, row_rank_x, 1);
        Uy = reshape(Uy, M, 1, row_rank_y);
        A.row_rank_ = row_ranges(4);
        A.col_rank_ = col_ranges(4);
        A.Umat_ = zeros(M, A.row_rank_);
        A.Vmat_ = zeros(A.col_size_, A.col_rank_);
        A.Umat_(:, 1 : row_ranges(2)) ...
            = reshape(Dx .* Uy, M, nx * row_rank_y);
        A.Umat_(:, (row_ranges(2) + 1) : row_ranges(3)) ...
            = reshape(Ux .* Dy, M, row_rank_x * ny);
        A.Umat_(:, (row_ranges(3) + 1) : row_ranges(4)) ...
            = reshape(Ux .* Uy, M, row_rank_x * row_rank_y);
        A.Vmat_(:, 1 : col_ranges(2)) = kron(Vy, eye(nx));
        A.Vmat_(:, (col_ranges(2) + 1) : col_ranges(3)) ...
            = kron(eye(ny), Vx);
        A.Vmat_(:, (col_ranges(3) + 1) : col_ranges(4)) ...
            = kron(Vy, Vx);
    end

    A.row_xy_ = [];
    return;
end

% Ranges of the three basis groups in each child.
child_row_ranges = zeros(A.num_children_, 4);
child_col_ranges = zeros(A.num_children_, 4);
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    if Ax.leaf_ == 1
        child_nx = Ax.col_size_;
        child_row_rank_x = Ax.row_rank_ * (Ax.level_ > 0);
        child_col_rank_x = Ax.col_rank_ * (Ax.level_ > 0);
    else
        child_nx = Ax.children_{ix}.col_size_;
        child_row_rank_x = Ax.children_{ix}.row_rank_;
        child_col_rank_x = Ax.children_{ix}.col_rank_;
    end
    if Ay.leaf_ == 1
        child_ny = Ay.col_size_;
        child_row_rank_y = Ay.row_rank_ * (Ay.level_ > 0);
        child_col_rank_y = Ay.col_rank_ * (Ay.level_ > 0);
    else
        child_ny = Ay.children_{iy}.col_size_;
        child_row_rank_y = Ay.children_{iy}.row_rank_;
        child_col_rank_y = Ay.children_{iy}.col_rank_;
    end
    child_row_ranges(i, :) = cumsum([0, ...
        child_nx * child_row_rank_y, ...
        child_row_rank_x * child_ny, ...
        child_row_rank_x * child_row_rank_y]);
    child_col_ranges(i, :) = cumsum([0, ...
        child_nx * child_col_rank_y, ...
        child_col_rank_x * child_ny, ...
        child_col_rank_x * child_col_rank_y]);
end

% Construct sibling interactions.
A.Bmat_ = cell(A.num_children_);
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    for j = [1 : (i - 1), (i + 1) : A.num_children_]
        jx = mod(j - 1, num_children_x) + 1;
        jy = floor((j - 1) / num_children_x) + 1;
        if ix == jx
            group = 1;
            K = kron(Ay.Bmat_{iy, jy}, ...
                eye(A.children_{i}.x_col_size_));
        elseif iy == jy
            group = 2;
            K = kron(eye(A.children_{i}.y_col_size_), ...
                Ax.Bmat_{ix, jx});
        else
            group = 3;
            K = kron(Ay.Bmat_{iy, jy}, Ax.Bmat_{ix, jx});
        end
        I = (child_row_ranges(i, group) + 1) ...
            : child_row_ranges(i, group + 1);
        J = (child_col_ranges(j, group) + 1) ...
            : child_col_ranges(j, group + 1);
        A.Bmat_{i, j} = zeros(A.children_{i}.row_rank_, ...
            A.children_{j}.col_rank_);
        A.Bmat_{i, j}(I, J) = K;
    end
end

if A.level_ == 0
    return;
end

% Construct transfer matrices.
A.row_rank_ = row_ranges(4);
A.col_rank_ = col_ranges(4);
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    if Ax.leaf_ == 1
        Ex = eye(nx);
        Rx = eye(row_rank_x);
        Wx = eye(col_rank_x);
        Fxi = zeros(row_rank_x, nx);
    else
        Ex = zeros(Ax.children_{ix}.col_size_, nx);
        J = Ax.children_{ix}.col_offset_ - Ax.col_offset_ ...
            + (1 : Ax.children_{ix}.col_size_);
        Ex(:, J) = eye(Ax.children_{ix}.col_size_);
        Rx = Ax.children_{ix}.Rmat_;
        Wx = Ax.children_{ix}.Wmat_;
        Fxi = Fx.mat{ix};
    end
    if Ay.leaf_ == 1
        Ey = eye(ny);
        Ry = eye(row_rank_y);
        Wy = eye(col_rank_y);
        Fyi = zeros(row_rank_y, ny);
    else
        Ey = zeros(Ay.children_{iy}.col_size_, ny);
        J = Ay.children_{iy}.col_offset_ - Ay.col_offset_ ...
            + (1 : Ay.children_{iy}.col_size_);
        Ey(:, J) = eye(Ay.children_{iy}.col_size_);
        Ry = Ay.children_{iy}.Rmat_;
        Wy = Ay.children_{iy}.Wmat_;
        Fyi = Fy.mat{iy};
    end

    A.children_{i}.Rmat_ = zeros(...
        A.children_{i}.row_rank_, A.row_rank_);
    A.children_{i}.Wmat_ = zeros(...
        A.children_{i}.col_rank_, A.col_rank_);

    I = 1 : child_row_ranges(i, 2);
    A.children_{i}.Rmat_(I, 1 : row_ranges(2)) = kron(Ry, Ex);
    I = (child_row_ranges(i, 2) + 1) : child_row_ranges(i, 3);
    A.children_{i}.Rmat_(I, (row_ranges(2) + 1) : row_ranges(3)) ...
        = kron(Ey, Rx);
    I = (child_row_ranges(i, 3) + 1) : child_row_ranges(i, 4);
    A.children_{i}.Rmat_(I, 1 : row_ranges(2)) = kron(Ry, Fxi);
    A.children_{i}.Rmat_(I, (row_ranges(2) + 1) : row_ranges(3)) ...
        = kron(Fyi, Rx);
    A.children_{i}.Rmat_(I, (row_ranges(3) + 1) : row_ranges(4)) ...
        = kron(Ry, Rx);

    J = 1 : child_col_ranges(i, 2);
    A.children_{i}.Wmat_(J, 1 : col_ranges(2)) = kron(Wy, Ex);
    J = (child_col_ranges(i, 2) + 1) : child_col_ranges(i, 3);
    A.children_{i}.Wmat_(J, (col_ranges(2) + 1) : col_ranges(3)) ...
        = kron(Ey, Wx);
    J = (child_col_ranges(i, 3) + 1) : child_col_ranges(i, 4);
    A.children_{i}.Wmat_(J, (col_ranges(3) + 1) : col_ranges(4)) ...
        = kron(Wy, Wx);
end

end
