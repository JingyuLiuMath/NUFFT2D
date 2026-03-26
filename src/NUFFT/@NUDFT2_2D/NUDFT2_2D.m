classdef NUDFT2_2D < handle
    
    % A[j, (kx, ky)] = exp(-2 * pi * 1i * x[j] * kx) ...
    %     * exp(-2 * pi * 1i * y[j] * ky),
    % 0 <= j < M, 0 <= kx < nx, 0 <= ky < ny.
    % N = nx * ny.
    
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