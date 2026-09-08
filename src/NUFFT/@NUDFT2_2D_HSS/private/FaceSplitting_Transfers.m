function [V, F] = FaceSplitting_Transfers(A)
% Cache F in D(parent)|child = D(child)*E + U(child)*F.

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

data.children = cell(1, A.num_children_);
V = zeros(A.col_size_, A.col_rank_ * (A.level_ > 0));
for i = 1 : A.num_children_
    [data.children{i}.V, F.children{i}] = FaceSplitting_Transfers(A.children_{i});
    if A.level_ > 0
        I = A.children_{i}.col_offset_ - A.col_offset_ + (1 : A.children_{i}.col_size_);
        V(I, :) = data.children{i}.V * A.children_{i}.Wmat_;
    end
end

if A.level_ > 0
    for i = 1 : A.num_children_
        F.mat{i} = zeros(A.children_{i}.row_rank_, A.col_size_);
        for j = 1 : A.num_children_
            if i ~= j
                J = A.children_{j}.col_offset_ - A.col_offset_ + (1 : A.children_{j}.col_size_);
                F.mat{i}(:, J) = A.Bmat_{i, j} * data.children{j}.V';
            end
        end
    end
end

end
