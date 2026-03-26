function f = URV_StarSolve(A, u)
% URV_StarSolve

% Jingyu Liu, December 9, 2024.

% Remark. In this case, we require that A be a square matrix.

arguments (Input)
    A HSS;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

A.FillVector_Col(u);

% Upward.
for level = A.max_level_ : -1 : 1
    A.URV_StarSolve_Upward(level);
    A.URV_StarSolve_Update(level);
end

% Root.
A.URV_StarSolve_Root();

% Downward
for level = 1 : 1 : A.max_level_
    A.URV_StarSolve_Downward(level);
end

f = A.FetchVector_Row();

end