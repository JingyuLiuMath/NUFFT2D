function u = URV_Solve(A, f)
% URV_Solve

% Jingyu Liu, November 22, 2024.

arguments (Input)
    A HSS;
    f (:, :) double;
end

arguments (Output)
    u (:, :) double;
end

A.FillVector_Row(f);

% Upward.
for level = A.max_level_ : -1 : 1
    A.URV_Solve_Upward(level);
end

% Root.
A.URV_Solve_Root();

% Downward
for level = 1 : 1 : A.max_level_
    A.URV_Solve_Downward(level);
end

u = A.FetchVector_Col();

end