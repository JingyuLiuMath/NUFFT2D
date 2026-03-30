classdef NUDFT2_HSS < HSS
    
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
        row_ind_ (:, 1) double;
        col_ind_ (:, 1) double;
        row_re_ (:, 1) double;
        col_re_ (:, 1) double;
        % -----------------------------------------------------------------
    end
    
     methods
         function A = NUDFT2_HSS(...
                 N, pos_start, pos_end, ...
                 level, row_offset, col_offset)
             arguments (Input)
                 N (1, 1) double;
                 pos_start (1, 1) double = 0;
                 pos_end (1, 1) double = N - 1;
                 level (1, 1) double = 0;
                 row_offset (1, 1) double = 0;
                 col_offset (1, 1) double = 0;
             end

             arguments (Output)
                 A NUDFT2_HSS;
             end

             A.col_global_size_ = N;
             A.pos_start_ = pos_start;
             A.pos_end_ = pos_end;
             A.col_size_ = pos_end - pos_start + 1;
             A.level_ = level;
             A.row_offset_ = row_offset;
             A.col_offset_ = col_offset;

         end
         
     end

end