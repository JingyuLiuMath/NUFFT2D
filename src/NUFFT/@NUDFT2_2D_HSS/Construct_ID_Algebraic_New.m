function Construct_ID_Algebraic_New(A, H, rank_or_tol)

arguments (Input)
    A NUDFT2_2D_HSS;
    H (:, :) double;
    rank_or_tol (1, 1) double;
end

for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_ID_Algebraic_New(...
        H, ...
        level, rank_or_tol);
    [row_act, col_act] = A.CollectActive(level);
    H = H(row_act, col_act);
end

A.ConstructRootGenerators4(H);

end