function mem = Storage(A)
% Storage

% Jingyu Liu, November 19, 2024.

arguments (Input)
    A HSS;
end

arguments (Output)
    mem (1, 1) double;
end

mem = numel(A.Amat_) + numel(A.Umat_) + numel(A.Vmat_);

for i = 1 : size(A.Rmat_, 2)
    mem = mem + numel(A.Rmat_{i});
end

for i = 1 : size(A.Wmat_, 2)
    mem = mem + numel(A.Wmat_{i});
end

for i = 1 : size(A.Bmat_, 1)
    for j = 1 : size(A.Bmat_, 2)
        mem = mem + numel(A.Bmat_{i, j});
    end
end

if A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        mem = mem + A.children_{i}.Storage();
    end
end

end