function BBC_FillAuxiliaryMatrix(A, Omega, Y, Psi, Z)
% BBC_FillAuxiliaryMatrix

% Jingyu Liu, December 5, 2024.

arguments (Input)
    A HSS;
    Omega (:, :) double;
    Y (:, :) double;
    Psi (:, :) double;
    Z (:, :) double;
end

if A.leaf_ == 1
    A.BBC_Omega_ = Omega;
    A.BBC_Y_ = Y;
    A.BBC_Psi_ = Psi;
    A.BBC_Z_ = Z;
else
    % Recursion.
    row_offset = 0;
    col_offset = 0;
    for i = 1 : A.num_children_
        current_row_size = A.children_{i}.row_size_;
        current_col_size = A.children_{i}.col_size_;
        A.children_{i}.BBC_FillAuxiliaryMatrix(...
            Omega((col_offset + 1) : (col_offset + current_col_size), :), ...
            Y((row_offset + 1) : (row_offset + current_row_size), :), ...
            Psi((row_offset + 1) : (row_offset + current_row_size), :), ...
            Z((col_offset + 1) : (col_offset + current_col_size), :));
        row_offset = row_offset + current_row_size;
        col_offset = col_offset + current_col_size;
    end
end

end