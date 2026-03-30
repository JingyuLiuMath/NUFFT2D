classdef NUDFT2_2D_HSS < HSS

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
        row_sk_size_ (1, 1) double;
        col_sk_size_ (1, 1) double;
        row_sk_offset_ (1, 1) double;
        col_sk_offset_(1, 1) double;
        row_ind_ (:, 1);
        col_ind_ (:, 1);
        row_re_ (:, 1);
        col_re_ (:, 1);
        % -----------------------------------------------------------------
    end
    
     methods
         function A = NUDFT2_2D_HSS(...
                 nx, ny, ...
                 x_pos_start, x_pos_end, ...
                 y_pos_start, y_pos_end, ...
                 level, row_offset, col_offset)
             arguments (Input)
                 nx (1, 1) double;
                 ny (1, 1) double
                 x_pos_start (1, 1) double = 0;
                 x_pos_end (1, 1) double = nx - 1;
                 y_pos_start (1, 1) double = 0;
                 y_pos_end (1, 1) double = ny - 1;
                 level (1, 1) double = 0;
                 row_offset (1, 1) double = 0;
                 col_offset (1, 1) double = 0;
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
             A.row_offset_ = row_offset;
             A.col_offset_ = col_offset;

         end
         
     end

end