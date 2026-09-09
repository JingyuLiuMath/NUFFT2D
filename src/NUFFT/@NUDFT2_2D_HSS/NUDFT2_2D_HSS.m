classdef NUDFT2_2D_HSS < HSS
    % NUDFT2_2D_HSS HSS representation of a transformed 2D NUDFT.
    % Construction data is local; all matrix properties are inherited.

    methods
        function A = NUDFT2_2D_HSS(nx, ny)
            arguments (Input)
                nx (1, 1) double;
                ny (1, 1) double;
            end

            A.col_global_size_ = nx * ny;
            A.col_size_ = nx * ny;
        end

        [p, q] = BuildTree(A, xy, nx, ny, N_leaf, ranges);
        Construct_ID_Proxy(A, xy, nx, ny, tol);
        [p, q, row_index] = BuildTree_FaceSplitting(A, Ax, Ay, px, py);
        Construct_FaceSplitting(A, Ax, Ay, row_index);
    end

end
