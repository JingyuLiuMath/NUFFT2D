function r = Rank(A)

arguments (Input)
    A HSS;
end

arguments (Output)
    r (1, 1) double;
end

r = max(A.row_rank_, A.col_rank_);

if A.leaf_ == 0
    for i = 1 : A.num_children_
        r = max(r, A.children_{i}.Rank());
    end
end

end
