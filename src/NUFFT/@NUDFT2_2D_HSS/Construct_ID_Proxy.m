function Construct_ID_Proxy(A, tol)
% Construct

% Jingyu Liu, February 4, 2025.

arguments (Input)
    A NUDFT2_2D_HSS;
    tol (1, 1) double;
end

for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_ID_Proxy(level, tol);
end

A.ConstructRootGenerators();

end