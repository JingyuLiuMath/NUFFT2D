function Construct_ID_Full(A, x, col_pos, rank_or_tol)

arguments (Input)
    A NUDFT2_HSS;
    x (:, 1) double;
    col_pos (:, 1) double;
    rank_or_tol (1, 1) double;
end

row_inact = [];
col_inact = [];
for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_ID_Full(...
        x, col_pos, ...
        row_inact, col_inact, ...
        level, rank_or_tol);
    [row_inact, col_inact] = A.Deactivate(level, row_inact, col_inact);
end

A.ConstructRootGenerators2(x, col_pos);

end