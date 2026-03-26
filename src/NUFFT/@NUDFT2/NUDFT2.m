classdef NUDFT2 < handle

    % A[j, k] = exp(-2 * pi * 1i * x[j] * k), 0 <= j < M, 0 <= k < N.

    % Row size: M.
    % Col size: N.
    
    properties
        M_ (1, 1) double;
        N_ (1, 1) double;
        AFinv_HSS_ NUDFT2_HSS;
        x_perm_ (:, 1) double;
        x_inv_perm_ (:, 1) double;
    end
    
     methods
         function A = NUDFT2(N)
             arguments (Input)
                 N (1, 1) double;
             end

             arguments (Output)
                 A NUDFT2;
             end

             A.N_ = N;
             
         end
         
     end

end