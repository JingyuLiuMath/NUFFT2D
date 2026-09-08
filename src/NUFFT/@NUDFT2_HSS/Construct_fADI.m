function Construct_fADI(A, x, tol)
% Construct_fADI Construct generators on the existing tree using x(p).
% Skeleton points are passed in a local struct and are not stored on A.

arguments (Input)
    A NUDFT2_HSS;
    x (:, 1) double;
    tol (1, 1) double;
end

ConstructGenerators_fADI(A, x, tol);

end
