function Construct_fADI(A, tol)
% Construct_fADI

% Jingyu Liu, April 14, 2025.

arguments (Input)
    A NUDFT2_HSS;
    tol (1, 1) double;
end

for level = A.max_level_ : -1 : 1
    A.ConstructGenerators_fADI(level, tol);
end

A.ConstructRootGenerators();

end