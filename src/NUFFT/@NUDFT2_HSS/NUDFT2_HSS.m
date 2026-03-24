classdef NUDFT2_HSS < HSS
    % NUDFT2_HSS

    % Jingyu Liu, November 18, 2024.

    % Let A be an M-by-N NUDFT2 matrx given by
    %     A(j, k) = exp(-2 * pi * 1i * x[j] * k), 0 <= j < M, 0 <= k < N.  
    % Suppose F is the N-by-N DFT matrix
    %     F(j, k) = exp(-2 * pi * 1i * j * k / N), 0 <= j, k < N.
    % Then tilde_A (= A * inv(F))'s entry has the representation formula:
    %    tilde_A(j, k) = k(z[j], w[k])
    % where k(z, w) = (z^N - 1) * w / (z - w) / N is the kernel function,
    % z[j] = gamma[j] = exp(-2 * pi * 1i * x[j]) and
    % w[k] = exp(-2 * pi * 1i * k / N).
    % tilde_A can be approximated by an HSS matrix.

    % Row size: M.
    % Col size: N.
    % Each row index j denotes a point exp(-2 * pi * 1i * x[j]).
    % Each col index k denotes a point exp(-2 * pi * 1i * k).

    % Note. We will always assume that N is the power of 2.
    
    properties
        % *****************************************************************
        % PROPERTY: Col Points.
        pos_start_ (1, 1) double;  % 0-base.
        pos_end_ (1, 1) double;  % 0-base.
        col_pos_ (:, 1) double;  % 0-base.
        % -----------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Row Points.
        row_x_ (:, 1) double;
        % -----------------------------------------------------------------
        
        % *****************************************************************
        % PROPERTY: HSS row and cols.
        row_ind_ (:, 1);
        col_ind_ (:, 1);
        self_row_ind_ (:, 1);
        self_col_ind_ (:, 1);
        % -----------------------------------------------------------------
    end
    
     methods
         function A = NUDFT2_HSS(...
                 N, pos_start, pos_end, ...
                 level, level_order, node_order)
             arguments (Input)
                 N (1, 1) double;
                 pos_start (1, 1) double = 0;
                 pos_end (1, 1) double = N - 1;
                 level (1, 1) double = 0;
                 level_order (1, 1) double = 0;
                 node_order (1, 1) double = 0;
             end

             arguments (Output)
                 A NUDFT2_HSS;
             end

             A.col_global_size_ = N;
             A.pos_start_ = pos_start;
             A.pos_end_ = pos_end;
             A.col_size_ = pos_end - pos_start + 1;
             A.level_ = level;
             A.level_order_ = level_order;
             A.node_order_ = node_order;

         end
         
     end

end