function Construct_ID_Proxy(A, xy, nx, ny, tol)
% Construct_ID_Proxy Construct generators on the existing tree using xy(p, :).

arguments (Input)
    A NUDFT2_2D_HSS;
    xy (:, 2) double;
    nx (1, 1) double;
    ny (1, 1) double;
    tol (1, 1) double;
end

ConstructGenerators_ID_Proxy(A, xy, nx, ny, [0, nx - 1, 0, ny - 1], tol);

end
