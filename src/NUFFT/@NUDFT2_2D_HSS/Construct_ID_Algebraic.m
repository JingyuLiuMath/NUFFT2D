function Construct_ID_Algebraic(A, H, rank_or_tol)

arguments (Input)
    A NUDFT2_2D_HSS;
    H(:, :) double;
    rank_or_tol (1, 1) double;
end

row_inact = [];
col_inact = [];
for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_ID_Algebraic(...
        H, ...
        row_inact, col_inact, ...
        level, rank_or_tol);
    [row_inact, col_inact] = A.Deactivate(level, row_inact, col_inact);
end

A.ConstructRootGenerators3(H);

end