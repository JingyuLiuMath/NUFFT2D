function Construct(A, rank_or_tol)
% Construct

% Jingyu Liu, February 4, 2025.

arguments (Input)
    A NUDFT2_2D_HSS;
    rank_or_tol (1, 1) double;
end

for level = A.max_level_ : -1 : 1
    A.ConstructGenerators(level, rank_or_tol);
end

A.ConstructRootGenerators();

end