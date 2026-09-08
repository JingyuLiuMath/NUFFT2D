classdef NUDFT2_HSS < HSS
    % NUDFT2_HSS HSS representation of a transformed 1D NUDFT.
    % The contiguous column interval is given by col_offset_ and col_size_.

    methods
        function A = NUDFT2_HSS(N)
            arguments (Input)
                N (1, 1) double;
            end

            A.col_global_size_ = N;
            A.col_size_ = N;
        end

        p = BuildTree(A, x, n_leaf);
        Construct_fADI(A, x, tol);
    end

end
