function Construct_ID_Proxy(A, rank_or_tol)
% Construct

% Jingyu Liu, November 24, 2024.

arguments (Input)
    A NUDFT2_HSS;
    rank_or_tol (1, 1) double;
end

for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_ID_Proxy(level, rank_or_tol);
end

A.ConstructRootGenerators();

end