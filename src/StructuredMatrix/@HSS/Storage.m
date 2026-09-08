function mem = Storage(A)
% Storage Total bytes of numeric generators in this subtree.

arguments (Input)
    A HSS;
end

arguments (Output)
    mem (1, 1) double;
end

mem = byte_size(A.Amat_) + byte_size(A.Umat_) + byte_size(A.Vmat_) ...
    + byte_size(A.Rmat_) + byte_size(A.Wmat_);

for i = 1 : size(A.Bmat_, 1)
    for j = 1 : size(A.Bmat_, 2)
        mem = mem + byte_size(A.Bmat_{i, j});
    end
end

if A.leaf_ == 0
    for i = 1 : A.num_children_
        mem = mem + A.children_{i}.Storage();
    end
end

end
