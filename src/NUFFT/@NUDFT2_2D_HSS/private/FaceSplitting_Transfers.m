function [F, V] = FaceSplitting_Transfers(A)
% D(parent)|child = D(child)*E + U(child)*F.mat{i}.

F.mat = cell(1, A.num_children_);
F.children = cell(1, A.num_children_);
if A.leaf_ == 1
    if A.level_ == 0
        V = zeros(A.col_size_, 0);
    else
        V = A.Vmat_;
    end
    return;
end

V_children = cell(1, A.num_children_);
V = zeros(A.col_size_, A.col_rank_ * (A.level_ > 0));
for i = 1 : A.num_children_
    [F.children{i}, V_children{i}] = FaceSplitting_Transfers(A.children_{i});
end
if A.level_ == 0
    return;
end

for i = 1 : A.num_children_
    I = A.children_{i}.col_offset_ - A.col_offset_ + (1 : A.children_{i}.col_size_);
    V(I, :) = V_children{i} * A.children_{i}.Wmat_;

    F.mat{i} = zeros(A.children_{i}.row_rank_, A.col_size_);
    for j = 1 : A.num_children_
        if i ~= j
            J = A.children_{j}.col_offset_ - A.col_offset_ + (1 : A.children_{j}.col_size_);
            F.mat{i}(:, J) = A.Bmat_{i, j} * V_children{j}';
        end
    end
end

end
