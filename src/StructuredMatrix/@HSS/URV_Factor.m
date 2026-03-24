function URV_Factor(A)
% URV

% Jingyu Liu, November 22, 2024.

arguments (Input)
    A HSS;
end

% Elimination and merge.
for level = A.max_level_ : -1 : 1
    A.URV_Eliminate(level);
end

A.URV_RootFactor();

end