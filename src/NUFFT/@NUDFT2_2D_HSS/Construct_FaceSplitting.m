function [p, q] = Construct_FaceSplitting(A, Ax, Ay, px, py)

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

% Build the product tree.
[p, q] = A.BuildTree_FaceSplitting(Ax, Ay, px, py);

% Cache only the additional matrices in D(parent)|child = D(child)*E + U(child)*F.
[~, Fx] = prepare_transfers(Ax);
[~, Fy] = prepare_transfers(Ay);

for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_FaceSplitting(level, Ax, Ay, px, py, Fx, Fy);
end

% Root.
A.ConstructRootGenerators_FaceSplitting(Ax, Ay, px, py);

end

function [V, F] = prepare_transfers(H)
F.mat = cell(1, H.num_children_);
F.children = cell(1, H.num_children_);
if H.leaf_
    V = H.Vmat_;
    if H.level_ == 0
        V = zeros(H.col_size_, 0);
    end
    return;
end

% Expanded column bases are used here and are not kept in the cache.
Vc = cell(1, H.num_children_);
V = zeros(H.col_size_, H.col_rank_ * (H.level_ > 0));
for i = 1 : H.num_children_
    C = H.children_{i};
    [Vc{i}, F.children{i}] = prepare_transfers(C);
    if H.level_ > 0
        I = C.col_offset_ - H.col_offset_ + (1 : C.col_size_);
        V(I, :) = Vc{i} * H.Wmat_{i};
    end
end

% The product root has no transfer matrices.
if H.level_ > 0
    for i = 1 : H.num_children_
        C = H.children_{i};
        F.mat{i} = zeros(C.row_rank_, H.col_size_);
        for j = [1 : (i - 1), (i + 1) : H.num_children_]
            D = H.children_{j};
            J = D.col_offset_ - H.col_offset_ + (1 : D.col_size_);
            F.mat{i}(:, J) = H.Bmat_{i, j} * Vc{j}';
        end
    end
end
end
