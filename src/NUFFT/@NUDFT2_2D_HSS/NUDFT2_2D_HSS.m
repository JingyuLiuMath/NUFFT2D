classdef NUDFT2_2D_HSS < HSS
    % NUDFT2_2D_HSS

    % Jingyu Liu, January 27, 2025.

    properties
        % *****************************************************************
        % PROPERTY: Col Points.
        nx_ (1, 1) double;
        ny_ (1, 1) double;
        x_pos_start_ (1, 1) double;  % 0-base.
        x_pos_end_ (1, 1) double;  % 0-base.
        x_col_pos_ (:, 1) double;  % 0-base.
        x_col_size_ (1, 1) double;
        y_pos_start_ (1, 1) double;  % 0-base.
        y_pos_end_ (1, 1) double;  % 0-base.
        y_col_pos_ (:, 1) double;  % 0-base.
        y_col_size_ (1, 1) double;
        col_pos_ (:, 2) double;  % 0-base.
        % -----------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Row Points.
        row_xy_ (:, 2) double;
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
         function A = NUDFT2_2D_HSS(...
                 nx, ny, ...
                 x_pos_start, x_pos_end, ...
                 y_pos_start, y_pos_end, ...
                 level, level_order, node_order)
             arguments (Input)
                 nx (1, 1) double;
                 ny (1, 1) double
                 x_pos_start (1, 1) double = 0;
                 x_pos_end (1, 1) double = nx - 1;
                 y_pos_start (1, 1) double = 0;
                 y_pos_end (1, 1) double = ny - 1;
                 level (1, 1) double = 0;
                 level_order (1, 1) double = 0;
                 node_order (1, 1) double = 0;
             end

             arguments (Output)
                 A NUDFT2_2D_HSS;
             end

             N = nx * ny;
             A.col_global_size_ = N;
             A.nx_ = nx;
             A.ny_ = ny;
             A.x_pos_start_ = x_pos_start;
             A.x_pos_end_ = x_pos_end;
             A.x_col_size_ = x_pos_end - x_pos_start + 1;
             A.y_pos_start_ = y_pos_start;
             A.y_pos_end_ = y_pos_end;
             A.y_col_size_ = y_pos_end - y_pos_start + 1;
             A.col_size_ = A.x_col_size_ * A.y_col_size_;
             A.level_ = level;
             A.level_order_ = level_order;
             A.node_order_ = node_order;

         end
         
     end

end