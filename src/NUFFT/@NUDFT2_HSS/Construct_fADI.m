function Construct_fADI(A, tol)
% Construct_fADI Construct generators from points stored by BuildTree.

arguments (Input)
    A NUDFT2_HSS;
    tol (1, 1) double;
end

ConstructGenerators_fADI(A, tol);

end
