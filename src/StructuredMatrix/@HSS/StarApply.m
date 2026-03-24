function u = StarApply(A, f)
% StarApply u = A' * f.

% Jingyu Liu, December 8, 2024.

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
    A.StarApply_Upward(level);
end

% Root.
A.StarApply_Root();

% Downward
for level = 1 : 1 : A.max_level_
    A.StarApply_Downward(level);
end

u = A.FetchVector_Col();

end