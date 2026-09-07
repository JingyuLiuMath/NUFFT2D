function [p, q] = BuildTree_FaceSplitting(A, Ax, Ay, px, py)
% Build the product tree directly from the two 1D leaf partitions.

arguments (Input)
    A NUDFT2_2D_HSS;
    Ax NUDFT2_HSS;
    Ay NUDFT2_HSS;
    px (:, 1) double;
    py (:, 1) double;
end

arguments (Output)
    p (:, 1) double;
    q (:, 1) double;
end

A.row_global_size_ = numel(px);
A.max_level_ = A.level_;
if Ax.leaf_ && Ay.leaf_
    A.leaf_ = 1;
    row_x = px(Ax.row_offset_ + (1 : Ax.row_size_));
    row_y = py(Ay.row_offset_ + (1 : Ay.row_size_));
    p = intersect(row_x, row_y, "sorted");
    q = reshape((Ax.pos_start_ : Ax.pos_end_)' ...
        + A.nx_ * (Ay.pos_start_ : Ay.pos_end_) + 1, [], 1);
else
    % Carry a finished axis unchanged while the other axis is split.
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
    A.num_children_ = numel(ix);
    A.children_ = cell(1, A.num_children_);
    p = cell(1, A.num_children_);
    q = p;
    row_offset = A.row_offset_;
    col_offset = A.col_offset_;
    for i = 1 : A.num_children_
        Cx = xc{ix(i)};
        Cy = yc{iy(i)};
        C = NUDFT2_2D_HSS(A.nx_, A.ny_, ...
            Cx.pos_start_, Cx.pos_end_, Cy.pos_start_, Cy.pos_end_, ...
            A.level_ + 1, row_offset, col_offset);
        [p{i}, q{i}] = C.BuildTree_FaceSplitting(Cx, Cy, px, py);
        A.children_{i} = C;
        row_offset = row_offset + C.row_size_;
        col_offset = col_offset + C.col_size_;
        A.max_level_ = max(A.max_level_, C.max_level_);
    end
    p = vertcat(p{:});
    q = vertcat(q{:});
end
A.row_size_ = numel(p);

end
