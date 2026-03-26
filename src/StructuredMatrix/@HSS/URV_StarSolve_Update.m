function URV_StarSolve_Update(A, max_level)
% URV_StarSolve_Update

% Jingyu Liu, December 10, 2024.

arguments (Input)
    A HSS;
    max_level (1, 1) double;
end

% Upward.
for level = max_level : -1 : 1
    A.StarSolve_Update_Upward(level, max_level);
end

% Root.
A.StarSolve_Update_Root();

% Downward
for level = 1 : 1 : max_level
    A.StarSolve_Update_Downward(level, max_level);
end

end