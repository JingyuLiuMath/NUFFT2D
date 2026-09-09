function data = ConstructGenerators_FaceSplitting(A, Ax, Ay, Fx, Fy, row_index)
% Construct the 2D HSS generators.

nx = Ax.col_size_;
ny = Ay.col_size_;
rx = Ax.row_rank_ * (Ax.level_ > 0);
ry = Ay.row_rank_ * (Ay.level_ > 0);
vx = Ax.col_rank_ * (Ax.level_ > 0);
vy = Ay.col_rank_ * (Ay.level_ > 0);
data.nx = nx;
data.ny = ny;
% Basis groups: y, x, xy.
data.row_ranges = cumsum([0, nx * ry, rx * ny, rx * ry]);
data.col_ranges = cumsum([0, nx * vy, vx * ny, vx * vy]);
A.row_rank_ = data.row_ranges(4);
A.col_rank_ = data.col_ranges(4);
I_y = 1 : data.row_ranges(2);
I_x = (data.row_ranges(2) + 1) : data.row_ranges(3);
I_xy = (data.row_ranges(3) + 1) : data.row_ranges(4);
J_y = 1 : data.col_ranges(2);
J_x = (data.col_ranges(2) + 1) : data.col_ranges(3);
J_xy = (data.col_ranges(3) + 1) : data.col_ranges(4);

% Leaf generators.
if A.leaf_ == 1
    M = A.row_size_;
    I = A.row_offset_ + (1 : M);
    rows_x = row_index.x(I) - Ax.row_offset_;
    rows_y = row_index.y(I) - Ay.row_offset_;
    Dx = reshape(Ax.Amat_(rows_x, :), M, nx, 1);
    Dy = reshape(Ay.Amat_(rows_y, :), M, 1, ny);
    A.Amat_ = reshape(Dx .* Dy, M, nx * ny);
    if A.level_ == 0
        return;
    end

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
    Ux = reshape(Ux, M, rx, 1);
    Uy = reshape(Uy, M, 1, ry);
    U = zeros(M, A.row_rank_);
    V = zeros(nx * ny, A.col_rank_);
    U(:, I_y) = reshape(Dx .* Uy, M, nx * ry);
    U(:, I_x) = reshape(Ux .* Dy, M, rx * ny);
    U(:, I_xy) = reshape(Ux .* Uy, M, rx * ry);
    V(:, J_y) = kron(Vy, eye(nx));
    V(:, J_x) = kron(eye(ny), Vx);
    V(:, J_xy) = kron(Vy, Vx);
    A.Umat_ = U;
    A.Vmat_ = V;
    return;
end

% Recursion.
num_children_x = max(Ax.num_children_, 1);
child_data = cell(1, A.num_children_);
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    if Ax.leaf_ == 1
        child_data{i} = ConstructGenerators_FaceSplitting( ...
            A.children_{i}, Ax, Ay.children_{iy}, Fx, Fy.children{iy}, row_index);
    elseif Ay.leaf_ == 1
        child_data{i} = ConstructGenerators_FaceSplitting( ...
            A.children_{i}, Ax.children_{ix}, Ay, Fx.children{ix}, Fy, row_index);
    else
        child_data{i} = ConstructGenerators_FaceSplitting( ...
            A.children_{i}, Ax.children_{ix}, Ay.children_{iy}, ...
            Fx.children{ix}, Fy.children{iy}, row_index);
    end
end

% Off-diagonal blocks. Each block uses one basis group.
A.Bmat_ = cell(A.num_children_);
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    for j = 1 : A.num_children_
        if i == j
            continue;
        end
        jx = mod(j - 1, num_children_x) + 1;
        jy = floor((j - 1) / num_children_x) + 1;
        if ix == jx
            group = 1;  % Only y changes.
            K = kron(Ay.Bmat_{iy, jy}, eye(child_data{i}.nx));
        elseif iy == jy
            group = 2;  % Only x changes.
            K = kron(eye(child_data{i}.ny), Ax.Bmat_{ix, jx});
        else
            group = 3;  % Both x and y change.
            K = kron(Ay.Bmat_{iy, jy}, Ax.Bmat_{ix, jx});
        end
        I = (child_data{i}.row_ranges(group) + 1) : child_data{i}.row_ranges(group + 1);
        J = (child_data{j}.col_ranges(group) + 1) : child_data{j}.col_ranges(group + 1);
        A.Bmat_{i, j} = zeros(A.children_{i}.row_rank_, A.children_{j}.col_rank_);
        A.Bmat_{i, j}(I, J) = K;
    end
end
if A.level_ == 0
    return;
end

% Transfer matrices.
for i = 1 : A.num_children_
    ix = mod(i - 1, num_children_x) + 1;
    iy = floor((i - 1) / num_children_x) + 1;
    % E selects child columns; F accounts for the other sibling columns.
    if Ax.leaf_ == 1
        Ex = eye(nx);
        Rx = eye(rx);
        Wx = eye(vx);
        Fxi = zeros(rx, nx);
    else
        Ex = zeros(Ax.children_{ix}.col_size_, nx);
        J = Ax.children_{ix}.col_offset_ - Ax.col_offset_ + (1 : Ax.children_{ix}.col_size_);
        Ex(:, J) = eye(Ax.children_{ix}.col_size_);
        Rx = Ax.children_{ix}.Rmat_;
        Wx = Ax.children_{ix}.Wmat_;
        Fxi = Fx.mat{ix};
    end
    if Ay.leaf_ == 1
        Ey = eye(ny);
        Ry = eye(ry);
        Wy = eye(vy);
        Fyi = zeros(ry, ny);
    else
        Ey = zeros(Ay.children_{iy}.col_size_, ny);
        J = Ay.children_{iy}.col_offset_ - Ay.col_offset_ + (1 : Ay.children_{iy}.col_size_);
        Ey(:, J) = eye(Ay.children_{iy}.col_size_);
        Ry = Ay.children_{iy}.Rmat_;
        Wy = Ay.children_{iy}.Wmat_;
        Fyi = Fy.mat{iy};
    end

    R = zeros(A.children_{i}.row_rank_, A.row_rank_);
    W = zeros(A.children_{i}.col_rank_, A.col_rank_);
    I = 1 : child_data{i}.row_ranges(2);
    R(I, I_y) = kron(Ry, Ex);
    I = (child_data{i}.row_ranges(2) + 1) : child_data{i}.row_ranges(3);
    R(I, I_x) = kron(Ey, Rx);
    I = (child_data{i}.row_ranges(3) + 1) : child_data{i}.row_ranges(4);
    R(I, I_y) = kron(Ry, Fxi);
    R(I, I_x) = kron(Fyi, Rx);
    R(I, I_xy) = kron(Ry, Rx);
    J = 1 : child_data{i}.col_ranges(2);
    W(J, J_y) = kron(Wy, Ex);
    J = (child_data{i}.col_ranges(2) + 1) : child_data{i}.col_ranges(3);
    W(J, J_x) = kron(Ey, Wx);
    J = (child_data{i}.col_ranges(3) + 1) : child_data{i}.col_ranges(4);
    W(J, J_xy) = kron(Wy, Wx);
    A.children_{i}.Rmat_ = R;
    A.children_{i}.Wmat_ = W;
end

end
