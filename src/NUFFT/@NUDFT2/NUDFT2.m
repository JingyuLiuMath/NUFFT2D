classdef NUDFT2 < handle
    % NUDFT2

    % Jingyu Liu, December 8, 2024.

    % NUDFT2 Problem.
    % f[j] = \sum_{k = 0}^{N - 1} exp(-2 * pi * 1i * x[j] * omega[k]),
    % j = 0, 1, ..., M - 1.
    % where x[j] are non-uniform and omega[k] = k + N_offset.

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