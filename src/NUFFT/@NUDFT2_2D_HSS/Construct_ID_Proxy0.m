function Construct_ID_Proxy0(A, rank_or_tol, xy, col_pos)
% Construct

% Jingyu Liu, February 4, 2025.

arguments (Input)
    A NUDFT2_2D_HSS;
    rank_or_tol (1, 1) double;
    xy (:, 2) double;
    col_pos (:, 2) double;
end

row_inact = [];
col_inact = [];
for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_ID_Proxy0(level, rank_or_tol, ...
        xy, col_pos, ...
        row_inact, col_inact);
    [row_inact, col_inact] = A.Deactivate(level, row_inact, col_inact);
end

A.ConstructRootGenerators2(xy, col_pos);

end