classdef NUDFT2_2D < handle
    % NUDFT2_2D

    % Jingyu Liu, January 27, 2025.

    % NUDFT2-2D Problem.
    % f(j) = \sum_{k1 = 0}^{N1 - 1} \sum_{k2 = 0}^{N2 - 1}
    %          exp(-2 * pi * 1i * (x(j, 1) * k1 + x(j, 2) *k2),
    % j = 0, 1, ..., M - 1.
    % where x[j] is non-uniform.

    % Row size: M.
    % Col size: N.
    
    properties
        M_ (1, 1) double;
        nx_ (1, 1) double;
        ny_ (1, 1) double;
        N_ (1, 1) double;
        AFinv_HSS_  NUDFT2_2D_HSS;
        xy_perm_ (:, 1) double;
        xy_inv_perm_ (:, 1) double;
        omega_perm_ (:, 1) double;
        omega_inv_perm_ (:, 1) double;
    end
    
     methods
         function A = NUDFT2_2D(nx, ny)
             arguments (Input)
                 nx (1, 1) double;
                 ny (1, 1) double;
             end

             arguments (Output)
                 A NUDFT2_2D;
             end

             A.nx_ = nx;
             A.ny_ = ny;
             A.N_ = nx * ny;
             
         end
         
     end

end