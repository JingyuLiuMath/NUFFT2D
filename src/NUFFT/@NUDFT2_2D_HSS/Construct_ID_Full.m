function Construct_ID_Full(A, xy, col_pos, rank_or_tol)
% Construct

% Jingyu Liu, November 24, 2024.

arguments (Input)
    A NUDFT2_2D_HSS;
    xy (:, 2) double;
    col_pos (:, 2) double;
    rank_or_tol (1, 1) double;
end

row_active = (1 : A.row_global_size_)';
col_active = (1 : A.col_global_size_)';
for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_ID_Full(...
        xy, col_pos, ...
        row_active, col_active, ...
        level, rank_or_tol);
    [row_active, col_active] = A.CollectActive(level);
end

A.ConstructRootGenerators2(xy, col_pos);

end