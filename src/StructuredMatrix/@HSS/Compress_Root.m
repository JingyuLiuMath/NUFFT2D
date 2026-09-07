function Compress_Root(A, tol)
% Update root couplings and compress the bases of its children.

arguments (Input)
    A HSS;
    tol (1, 1) double;
end

if A.leaf_ == 1
    return;
end

for i = 1 : A.num_children_
    for j = 1 : A.num_children_
        if i ~= j
            A.Bmat_{i, j} = A.children_{i}.Compress_Tu_ ...
                * A.Bmat_{i, j} * A.children_{j}.Compress_Tv_';
        else
            A.Bmat_{i, j} = [];
        end
    end
end
% Clear after all sibling couplings have used the child transforms.
for i = 1 : A.num_children_
    A.children_{i}.Compress_Tu_ = [];
    A.children_{i}.Compress_Tv_ = [];
end

P = cell(1, A.num_children_);
Q = P;
% Form all sibling bases before changing any coupling matrices.
for i = 1 : A.num_children_
    Bu = cell(1, A.num_children_);
    Bv = Bu;
    Bu{i} = zeros(A.children_{i}.row_rank_, 0);
    Bv{i} = zeros(A.children_{i}.col_rank_, 0);
    for j = 1 : A.num_children_
        if i ~= j
            Bu{j} = A.Bmat_{i, j};
            Bv{j} = A.Bmat_{j, i}';
        end
    end
    [P{i}, A.children_{i}.Compress_Su_] = MySVDSketch(horzcat(Bu{:}), tol);
    [Q{i}, A.children_{i}.Compress_Sv_] = MySVDSketch(horzcat(Bv{:}), tol);
end

for i = 1 : A.num_children_
    for j = 1 : A.num_children_
        if i ~= j
            A.Bmat_{i, j} = P{i}' * A.Bmat_{i, j} * Q{j};
        end
    end
    C = A.children_{i};
    C.row_rank_ = size(P{i}, 2);
    C.col_rank_ = size(Q{i}, 2);
    if C.leaf_ == 1
        C.Umat_ = C.Umat_ * P{i};
        C.Vmat_ = C.Vmat_ * Q{i};
    else
        % Pass the parent basis changes down, as in (5.10).
        for j = 1 : C.num_children_
            C.Rmat_{j} = C.Rmat_{j} * P{i};
            C.Wmat_{j} = C.Wmat_{j} * Q{i};
        end
    end
end

end
