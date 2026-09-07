function ConstructRootGenerators_FaceSplitting(A, Ax, Ay, px, py)

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    px (:, 1) double;
    py (:, 1) double;
end

if A.leaf_ == 1
    nx = Ax.col_size_;
    ny = Ay.col_size_;
    row_x = px(Ax.row_offset_ + (1 : Ax.row_size_));
    row_y = py(Ay.row_offset_ + (1 : Ay.row_size_));
    [~, row_x, row_y] = intersect(row_x, row_y, "sorted");
    M = A.row_size_;
    Dx = reshape(Ax.Amat_(row_x, :), M, nx, 1);
    Dy = reshape(Ay.Amat_(row_y, :), M, 1, ny);
    A.Amat_ = reshape(Dx .* Dy, M, nx * ny);
else
    xc = Ax.children_;
    yc = Ay.children_;
    if Ax.leaf_
        xc = {Ax};
    end
    if Ay.leaf_
        yc = {Ay};
    end
    [ix, iy] = ndgrid(1 : numel(xc), 1 : numel(yc));
    ix = ix(:);
    iy = iy(:);

    row_ranges = zeros(A.num_children_, 4);
    col_ranges = zeros(A.num_children_, 4);
    for i = 1 : A.num_children_
        Cx = xc{ix(i)};
        Cy = yc{iy(i)};
        rx = Cx.row_rank_ * (Cx.level_ > 0);
        ry = Cy.row_rank_ * (Cy.level_ > 0);
        vx = Cx.col_rank_ * (Cx.level_ > 0);
        vy = Cy.col_rank_ * (Cy.level_ > 0);
        row_ranges(i, :) = cumsum([0, Cx.col_size_ * ry, rx * Cy.col_size_, rx * ry]);
        col_ranges(i, :) = cumsum([0, Cx.col_size_ * vy, vx * Cy.col_size_, vx * vy]);
    end

    A.Bmat_ = cell(A.num_children_);
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
    end
    % Clear after all root couplings have used the child factors.
    for i = 1 : A.num_children_
        A.children_{i}.FaceSplitting_Tu_ = [];
        A.children_{i}.FaceSplitting_Tv_ = [];
    end
end

end
